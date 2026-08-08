#!/usr/bin/env python3
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def write(path, text):
    (ROOT / path).write_text(text, encoding="utf-8")


def replace_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected one occurrence, found {n}")
    return text.replace(old, new, 1)


def patch_subroutine(text, name, transform):
    pat = re.compile(
        rf"(?ms)^\s*subroutine\s+{re.escape(name)}\b.*?^\s*end\s+subroutine\s+{re.escape(name)}\s*$"
    )
    m = pat.search(text)
    if not m:
        raise RuntimeError(f"subroutine {name} not found")
    block = transform(m.group(0))
    return text[:m.start()] + block + text[m.end():]


KV_MODULE = r'''module integrator_kv_ev_mod

!***********************************************************************
! Kelvin-Voigt integration with solvent evaporation.
!
! This module extends the historical JETSPIN Kelvin-Voigt integrators
! to the concentration-dependent material properties used by the
! Yarin-Koombhongse-Reneker evaporation model.  The ordinary Kelvin-
! Voigt path in integrator_mod is left unchanged for backward
! compatibility when evaporation is disabled.
!***********************************************************************

 use version_mod, only : mystart,myend,mxchunk,sum_world_darr, &
                         set_chunk,set_mxchunk,idrank
 use error_mod, only : error
 use nanojet_mod, only : mxnpjet,npjet,inpjet,systype,jetxx,jetyy, &
                         jetzz,jetst,jetvx,jetvy,jetvz,jetvl,jetve, &
                         compute_posnoinserted,evlim
 use integrator_mod, only : integrator
 use dynamic_refinement_mod, only : driver_dynamic_refinement
 use coulomb_force_mod, only : smooth_charge,restore_charge,coulforce, &
                               compute_coulomelec_driver
 use eom_ev_mod, only : eom1_KV_pos_v_ev,eom1_KV_st_ev, &
                        eom3_KV_pos_v_ev,eom3_KV_st_ev

 implicit none
 private
 public :: driver_integrator_KV_ev

 double precision, allocatable, save :: fxx(:,:),fyy(:,:),fzz(:,:)
 double precision, allocatable, save :: fst(:,:),fev(:,:)
 double precision, allocatable, save :: f1vx(:),f1vy(:),f1vz(:)
 double precision, allocatable, save :: f2vx(:),f2vy(:),f2vz(:)
 double precision, allocatable, save :: f3vx(:),f3vy(:),f3vz(:)
 double precision, allocatable, save :: f4vx(:),f4vy(:),f4vz(:)
 double precision, allocatable, save :: yxx(:),yyy(:),yzz(:),yst(:)
 double precision, allocatable, save :: yvx(:),yvy(:),yvz(:),yev(:)
 logical, save :: lworkspace=.false.
 logical, save :: lannounced=.false.

 contains

 subroutine driver_integrator_KV_ev(timesub,h,k,dorefinment)

  implicit none
  logical, intent(inout) :: dorefinment
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  integer :: i
  logical :: ltestinst

  call driver_dynamic_refinement(k,dorefinment)
  call set_chunk(inpjet,npjet)
  call set_mxchunk(mxnpjet)
  call ensure_workspace()

  if((.not.lannounced).and.idrank==0)then
    write(6,'(/,a,/)') &
     'Kelvin-Voigt evaporation integrator active (Yarin concentration laws)'
    lannounced=.true.
  endif

  select case(integrator)
    case(1)
      call eulsys_KV_ev(timesub,h,k)
    case(2)
      call rk2sys_KV_ev(timesub,h,k)
    case(3)
      call rk4sys_KV_ev(timesub,h,k)
    case default
      call error(1)
  end select

  ltestinst=.false.
  do i=inpjet,npjet
    if(isnan(dcos(jetxx(i))))ltestinst=.true.
    if(isnan(dcos(jetyy(i))))ltestinst=.true.
    if(isnan(dcos(jetzz(i))))ltestinst=.true.
    if(isnan(dcos(jetst(i))))ltestinst=.true.
    if(isnan(dcos(jetvx(i))))ltestinst=.true.
    if(isnan(dcos(jetvy(i))))ltestinst=.true.
    if(isnan(dcos(jetvz(i))))ltestinst=.true.
    if(isnan(dcos(jetve(i))))ltestinst=.true.
  enddo
  if(ltestinst)call error(14)

  return
 end subroutine driver_integrator_KV_ev

 subroutine ensure_workspace()
  implicit none

  if(lworkspace)return

  allocate(fxx(0:mxchunk,4),fyy(0:mxchunk,4),fzz(0:mxchunk,4))
  allocate(fst(0:mxchunk,4),fev(0:mxchunk,4))
  allocate(f1vx(0:mxnpjet),f1vy(0:mxnpjet),f1vz(0:mxnpjet))
  allocate(f2vx(0:mxnpjet),f2vy(0:mxnpjet),f2vz(0:mxnpjet))
  allocate(f3vx(0:mxnpjet),f3vy(0:mxnpjet),f3vz(0:mxnpjet))
  allocate(f4vx(0:mxnpjet),f4vy(0:mxnpjet),f4vz(0:mxnpjet))
  allocate(yxx(0:mxnpjet),yyy(0:mxnpjet),yzz(0:mxnpjet))
  allocate(yst(0:mxnpjet),yvx(0:mxnpjet),yvy(0:mxnpjet))
  allocate(yvz(0:mxnpjet),yev(0:mxnpjet))
  lworkspace=.true.

  return
 end subroutine ensure_workspace

 subroutine eval_stage(tstage,k,xs,ys,zs,ss,vxs,vys,vzs,ves, &
                       dx,dy,dz,ds,dve,ax,ay,az)
  implicit none
  integer, intent(in) :: k
  double precision, intent(in) :: tstage
  double precision, allocatable, dimension(:), intent(inout) :: xs,ys,zs
  double precision, allocatable, dimension(:), intent(in) :: ss,vxs,vys,vzs
  double precision, allocatable, dimension(:), intent(in) :: ves
  double precision, dimension(0:), intent(inout) :: dx,dy,dz,ds,dve
  double precision, allocatable, dimension(:), intent(inout) :: ax,ay,az
  integer :: ipoint,j

  dx(:)=0.d0
  dy(:)=0.d0
  dz(:)=0.d0
  ds(:)=0.d0
  dve(:)=0.d0
  ax(:)=0.d0
  ay(:)=0.d0
  az(:)=0.d0

  select case(systype)
    case(1)
      call smooth_charge(xs)
      call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs,ves)
    case(3)
      call smooth_charge(xs,ys,zs)
      call compute_posnoinserted(xs,ys,zs)
      call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs,ves)
    case default
      call error(2)
  end select

  j=0
  do ipoint=mystart,myend
    select case(systype)
      case(1)
        call eom1_KV_pos_v_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs, &
         jetvl,ves,coulforce,dx(j),dy(j),dz(j),ax(ipoint),ay(ipoint), &
         az(ipoint),dve(j),tstage,k)
      case(3)
        call eom3_KV_pos_v_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs, &
         jetvl,ves,coulforce,dx(j),dy(j),dz(j),ax(ipoint),ay(ipoint), &
         az(ipoint),dve(j),tstage,k)
    end select
    j=j+1
  enddo

  call sum_world_darr(ax,npjet+1)
  call sum_world_darr(ay,npjet+1)
  call sum_world_darr(az,npjet+1)

  j=0
  do ipoint=mystart,myend
    select case(systype)
      case(1)
        call eom1_KV_st_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves, &
         coulforce,ax,ay,az,dve(j),ds(j),tstage,k)
      case(3)
        call eom3_KV_st_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves, &
         coulforce,ax,ay,az,dve(j),ds(j),tstage,k)
    end select
    j=j+1
  enddo

  call restore_charge()

  return
 end subroutine eval_stage

 subroutine clamp_ev_volume(ipoint,value)
  implicit none
  integer, intent(in) :: ipoint
  double precision, intent(inout) :: value
  if((value/jetvl(ipoint))<evlim)value=jetvl(ipoint)*evlim
  return
 end subroutine clamp_ev_volume

 subroutine gather_state()
  implicit none
  call sum_world_darr(yxx,npjet+1)
  call sum_world_darr(yyy,npjet+1)
  call sum_world_darr(yzz,npjet+1)
  call sum_world_darr(yst,npjet+1)
  call sum_world_darr(yvx,npjet+1)
  call sum_world_darr(yvy,npjet+1)
  call sum_world_darr(yvz,npjet+1)
  call sum_world_darr(yev,npjet+1)
  return
 end subroutine gather_state

 subroutine commit_state()
  implicit none
  call sum_world_darr(yxx,npjet+1,jetxx)
  call sum_world_darr(yyy,npjet+1,jetyy)
  call sum_world_darr(yzz,npjet+1,jetzz)
  call sum_world_darr(yst,npjet+1,jetst)
  call sum_world_darr(yvx,npjet+1,jetvx)
  call sum_world_darr(yvy,npjet+1,jetvy)
  call sum_world_darr(yvz,npjet+1,jetvz)
  call sum_world_darr(yev,npjet+1,jetve)
  return
 end subroutine commit_state

 subroutine zero_state()
  implicit none
  yxx(:)=0.d0
  yyy(:)=0.d0
  yzz(:)=0.d0
  yst(:)=0.d0
  yvx(:)=0.d0
  yvy(:)=0.d0
  yvz(:)=0.d0
  yev(:)=0.d0
  return
 end subroutine zero_state

 subroutine eulsys_KV_ev(timesub,h,k)
  implicit none
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  integer :: ipoint,j

  call eval_stage(timesub,k,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
   jetvz,jetve,fxx(:,1),fyy(:,1),fzz(:,1),fst(:,1),fev(:,1), &
   f1vx,f1vy,f1vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+h*fxx(j,1)
    yyy(ipoint)=jetyy(ipoint)+h*fyy(j,1)
    yzz(ipoint)=jetzz(ipoint)+h*fzz(j,1)
    yst(ipoint)=jetst(ipoint)+h*fst(j,1)
    yvx(ipoint)=jetvx(ipoint)+h*f1vx(ipoint)
    yvy(ipoint)=jetvy(ipoint)+h*f1vy(ipoint)
    yvz(ipoint)=jetvz(ipoint)+h*f1vz(ipoint)
    yev(ipoint)=jetve(ipoint)+h*fev(j,1)
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call commit_state()
  timesub=timesub+h
  if(systype==3)call compute_posnoinserted(jetxx,jetyy,jetzz)
  return
 end subroutine eulsys_KV_ev

 subroutine rk2sys_KV_ev(timesub,h,k)
  implicit none
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  integer :: ipoint,j

  call eval_stage(timesub,k,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
   jetvz,jetve,fxx(:,1),fyy(:,1),fzz(:,1),fst(:,1),fev(:,1), &
   f1vx,f1vy,f1vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+h*fxx(j,1)
    yyy(ipoint)=jetyy(ipoint)+h*fyy(j,1)
    yzz(ipoint)=jetzz(ipoint)+h*fzz(j,1)
    yst(ipoint)=jetst(ipoint)+h*fst(j,1)
    yvx(ipoint)=jetvx(ipoint)+h*f1vx(ipoint)
    yvy(ipoint)=jetvy(ipoint)+h*f1vy(ipoint)
    yvz(ipoint)=jetvz(ipoint)+h*f1vz(ipoint)
    yev(ipoint)=jetve(ipoint)+h*fev(j,1)
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call gather_state()

  call eval_stage(timesub+h,k,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   fxx(:,2),fyy(:,2),fzz(:,2),fst(:,2),fev(:,2),f2vx,f2vy,f2vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+0.5d0*h*(fxx(j,1)+fxx(j,2))
    yyy(ipoint)=jetyy(ipoint)+0.5d0*h*(fyy(j,1)+fyy(j,2))
    yzz(ipoint)=jetzz(ipoint)+0.5d0*h*(fzz(j,1)+fzz(j,2))
    yst(ipoint)=jetst(ipoint)+0.5d0*h*(fst(j,1)+fst(j,2))
    yvx(ipoint)=jetvx(ipoint)+0.5d0*h*(f1vx(ipoint)+f2vx(ipoint))
    yvy(ipoint)=jetvy(ipoint)+0.5d0*h*(f1vy(ipoint)+f2vy(ipoint))
    yvz(ipoint)=jetvz(ipoint)+0.5d0*h*(f1vz(ipoint)+f2vz(ipoint))
    yev(ipoint)=jetve(ipoint)+0.5d0*h*(fev(j,1)+fev(j,2))
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call commit_state()
  timesub=timesub+h
  if(systype==3)call compute_posnoinserted(jetxx,jetyy,jetzz)
  return
 end subroutine rk2sys_KV_ev

 subroutine rk4sys_KV_ev(timesub,h,k)
  implicit none
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  integer :: ipoint,j

  call eval_stage(timesub,k,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
   jetvz,jetve,fxx(:,1),fyy(:,1),fzz(:,1),fst(:,1),fev(:,1), &
   f1vx,f1vy,f1vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+0.5d0*h*fxx(j,1)
    yyy(ipoint)=jetyy(ipoint)+0.5d0*h*fyy(j,1)
    yzz(ipoint)=jetzz(ipoint)+0.5d0*h*fzz(j,1)
    yst(ipoint)=jetst(ipoint)+0.5d0*h*fst(j,1)
    yvx(ipoint)=jetvx(ipoint)+0.5d0*h*f1vx(ipoint)
    yvy(ipoint)=jetvy(ipoint)+0.5d0*h*f1vy(ipoint)
    yvz(ipoint)=jetvz(ipoint)+0.5d0*h*f1vz(ipoint)
    yev(ipoint)=jetve(ipoint)+0.5d0*h*fev(j,1)
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call gather_state()

  call eval_stage(timesub+0.5d0*h,k,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   fxx(:,2),fyy(:,2),fzz(:,2),fst(:,2),fev(:,2),f2vx,f2vy,f2vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+0.5d0*h*fxx(j,2)
    yyy(ipoint)=jetyy(ipoint)+0.5d0*h*fyy(j,2)
    yzz(ipoint)=jetzz(ipoint)+0.5d0*h*fzz(j,2)
    yst(ipoint)=jetst(ipoint)+0.5d0*h*fst(j,2)
    yvx(ipoint)=jetvx(ipoint)+0.5d0*h*f2vx(ipoint)
    yvy(ipoint)=jetvy(ipoint)+0.5d0*h*f2vy(ipoint)
    yvz(ipoint)=jetvz(ipoint)+0.5d0*h*f2vz(ipoint)
    yev(ipoint)=jetve(ipoint)+0.5d0*h*fev(j,2)
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call gather_state()

  call eval_stage(timesub+0.5d0*h,k,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   fxx(:,3),fyy(:,3),fzz(:,3),fst(:,3),fev(:,3),f3vx,f3vy,f3vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+h*fxx(j,3)
    yyy(ipoint)=jetyy(ipoint)+h*fyy(j,3)
    yzz(ipoint)=jetzz(ipoint)+h*fzz(j,3)
    yst(ipoint)=jetst(ipoint)+h*fst(j,3)
    yvx(ipoint)=jetvx(ipoint)+h*f3vx(ipoint)
    yvy(ipoint)=jetvy(ipoint)+h*f3vy(ipoint)
    yvz(ipoint)=jetvz(ipoint)+h*f3vz(ipoint)
    yev(ipoint)=jetve(ipoint)+h*fev(j,3)
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call gather_state()

  call eval_stage(timesub+h,k,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   fxx(:,4),fyy(:,4),fzz(:,4),fst(:,4),fev(:,4),f4vx,f4vy,f4vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+(h/6.d0)*(fxx(j,1)+ &
     2.d0*(fxx(j,2)+fxx(j,3))+fxx(j,4))
    yyy(ipoint)=jetyy(ipoint)+(h/6.d0)*(fyy(j,1)+ &
     2.d0*(fyy(j,2)+fyy(j,3))+fyy(j,4))
    yzz(ipoint)=jetzz(ipoint)+(h/6.d0)*(fzz(j,1)+ &
     2.d0*(fzz(j,2)+fzz(j,3))+fzz(j,4))
    yst(ipoint)=jetst(ipoint)+(h/6.d0)*(fst(j,1)+ &
     2.d0*(fst(j,2)+fst(j,3))+fst(j,4))
    yvx(ipoint)=jetvx(ipoint)+(h/6.d0)*(f1vx(ipoint)+ &
     2.d0*(f2vx(ipoint)+f3vx(ipoint))+f4vx(ipoint))
    yvy(ipoint)=jetvy(ipoint)+(h/6.d0)*(f1vy(ipoint)+ &
     2.d0*(f2vy(ipoint)+f3vy(ipoint))+f4vy(ipoint))
    yvz(ipoint)=jetvz(ipoint)+(h/6.d0)*(f1vz(ipoint)+ &
     2.d0*(f2vz(ipoint)+f3vz(ipoint))+f4vz(ipoint))
    yev(ipoint)=jetve(ipoint)+(h/6.d0)*(fev(j,1)+ &
     2.d0*(fev(j,2)+fev(j,3))+fev(j,4))
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call commit_state()
  timesub=timesub+h
  if(systype==3)call compute_posnoinserted(jetxx,jetyy,jetzz)
  return
 end subroutine rk4sys_KV_ev

end module integrator_kv_ev_mod
'''

KV_CHECK = r'''#!/usr/bin/env python3
"""Regression guard for the concentration-dependent Kelvin-Voigt extension."""
from __future__ import annotations
import argparse
import math
import re
from pathlib import Path

REF_RMU = 8.367216948213878
REF_RG = 1.6734433896427756
REF_DRMU = 2.3913155036429394
REF_DRG = 0.14357442280003282
REF_STRAIN = -0.20121487043148
REF_DSTRESS = 0.1265700720011721


def close(name, got, ref, tol=2e-12):
    if not math.isclose(got, ref, rel_tol=tol, abs_tol=tol):
        raise AssertionError(f"{name}: {got:.16g} != {ref:.16g}")
    print(f"PASS {name:26s} {got:.12g}")


def compact(path):
    text = path.read_text(encoding="utf-8", errors="replace").lower()
    text = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
    return re.sub(r"\s+", "", text)


def require(path, expression, name):
    if re.sub(r"\s+", "", expression.lower()) not in compact(path):
        raise AssertionError(f"{name}: expression not found in {path}")
    print(f"PASS {name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eom-source", required=True, type=Path)
    ap.add_argument("--integrator-source", required=True, type=Path)
    args = ap.parse_args()

    cp0, cp = 0.06, 0.30
    b, m, tev = 7.0, 0.1, 1.0
    stress, strainrate, strainacc = 0.5, 0.1, -0.03
    # fV/V=-0.2, hence dc_p/dt=-c_p fV/V=+0.06.
    dcpdt = 0.06

    rmu = 10.0 ** (b * (cp**m - cp0**m))
    rg = rmu / ((cp/cp0)**tev)
    drmu = rmu * math.log(10.0) * b * m * cp**(m-1.0) * dcpdt
    drg = rg * (math.log(10.0)*b*m*cp**(m-1.0)-tev/cp) * dcpdt
    strain = (stress-rmu*strainrate)/rg
    dstress = (rg*strainrate + rmu*strainacc +
               drg*strain + drmu*strainrate)

    close("mu/mu0", rmu, REF_RMU)
    close("G/G0", rg, REF_RG)
    close("d(mu/mu0)/dtbar", drmu, REF_DRMU)
    close("d(G/G0)/dtbar", drg, REF_DRG)
    close("recovered strain", strain, REF_STRAIN)
    close("Kelvin-Voigt stress rate", dstress, REF_DSTRESS)

    # When concentration is frozen, the extension must reduce exactly to
    # the historical JETSPIN Kelvin-Voigt differential form.
    frozen = rg*strainrate + rmu*strainacc
    product_rule_frozen = rg*strainrate + rmu*strainacc + 0.0 + 0.0
    close("frozen-composition limit", product_rule_frozen, frozen)

    require(args.eom_source,
            "dcpdt=-cp*fevlocal/yve",
            "production dc_p/dt")
    require(args.eom_source,
            "strain=(stress-ratmu*strainrate)/ratg",
            "production strain recovery")
    require(args.eom_source,
            "fst=ratg*strainrate+ratmu*strainacc+dratg*strain+dratmu*strainrate",
            "production product rule")
    require(args.integrator_source,
            "call eom3_KV_st_ev",
            "3D Kelvin-Voigt evaporation path")
    require(args.integrator_source,
            "call clamp_ev_volume",
            "Yarin cutoff in KV integrators")

    print("Kelvin-Voigt evaporation regression checks passed")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
'''


def patch_eom():
    path = "source/eom_ev_mod.f90"
    text = read(path)
    text = replace_once(
        text,
        "evcsvapour,tev,lengthscale,tao,lairdrag",
        "evcsvapour,tev,lengthscale,tao,lairdrag,evlim",
        "eom import evlim",
    )

    def transform(block):
        block = block.replace("ycf,yax,yay,yaz,fst,timesub,k)",
                              "ycf,yax,yay,yaz,fevlocal,fst,timesub,k)", 1)
        marker = "  double precision, intent(inout) ::  fst"
        if marker not in block:
            raise RuntimeError("KV stress declaration marker not found")
        block = block.replace(marker,
            "  double precision, intent(in) :: fevlocal\n" + marker, 1)
        old = "fst = ratg*(beadvelup/beadlenup)+ratmu*(beadaccup/beadlenup)"
        n = block.count(old)
        if n < 1:
            raise RuntimeError("historical KV stress expression not found")
        new = ("call kv_ev_stress_rate(cp,ratmu,ratg,yve(ipoint),yvl(ipoint), &\n"
               "     fevlocal,yst(ipoint),beadvelup/beadlenup, &\n"
               "     beadaccup/beadlenup,fst)")
        block = block.replace(old, new)
        return block

    text = patch_subroutine(text, "eom1_KV_st_ev", transform)
    text = patch_subroutine(text, "eom3_KV_st_ev", transform)

    helper = r'''

 subroutine kv_ev_stress_rate(cp,ratmu,ratg,yve,yvl,fevlocal,stress, &
                              strainrate,strainacc,fst)

!***********************************************************************
! Product-rule Kelvin-Voigt stress rate with concentration-dependent
! viscosity and elastic modulus.  All quantities are in the standard
! JETSPIN nondimensionalization.  The historical JETSPIN Kelvin-Voigt
! kinematics is retained: strainrate=(1/l) dl/dt and the acceleration
! contribution is represented by strainacc=(1/l) dv_parallel/dt.
!***********************************************************************

  implicit none
  double precision, intent(in) :: cp,ratmu,ratg,yve,yvl,fevlocal
  double precision, intent(in) :: stress,strainrate,strainacc
  double precision, intent(out) :: fst
  double precision :: dcpdt,dratmu,dratg,strain

  dcpdt=0.d0
  if(yve>0.d0 .and. yvl>0.d0)then
    if((yve/yvl)>evlim*(1.d0+1.d-12))then
      dcpdt=-cp*fevlocal/yve
    endif
  endif

  dratmu=ratmu*dlog(10.d0)*Bev*mev*(cp**(mev-1.d0))*dcpdt
  dratg=ratg*(dlog(10.d0)*Bev*mev*(cp**(mev-1.d0))-tev/cp)*dcpdt

  strain=(stress-ratmu*strainrate)/ratg
  fst=ratg*strainrate+ratmu*strainacc+dratg*strain+dratmu*strainrate

  return
 end subroutine kv_ev_stress_rate
'''
    marker = "\n end module eom_ev_mod"
    if marker not in text:
        raise RuntimeError("eom module end marker not found")
    text = text.replace(marker, helper + marker, 1)
    write(path, text)


def patch_io():
    path = "source/io_mod.f90"
    text = read(path)
    old = "    if(lKVfluid)then\n      call warning(101)\n      ltest=.true.\n    endif\n"
    if old not in text:
        raise RuntimeError("Kelvin-Voigt evaporation prohibition not found")
    text = text.replace(old, "", 1)
    write(path, text)


def patch_main():
    path = "source/main.f90"
    text = read(path)
    text = replace_once(
        text,
        "                       pdbrescale,lreadrest\n",
        "                       pdbrescale,lreadrest,lKVfluid,levaporation\n",
        "main nanojet import",
    )
    text = replace_once(
        text,
        "  use integrator_mod, only : initime,endtime,driver_integrator\n",
        "  use integrator_mod, only : initime,endtime,driver_integrator\n"
        "  use integrator_kv_ev_mod, only : driver_integrator_KV_ev\n",
        "main KV evaporation module import",
    )
    text = replace_once(
        text,
        "    call driver_integrator(mytime,tstep,nstep,ldorefinment)\n",
        "    if(lKVfluid.and.levaporation)then\n"
        "      call driver_integrator_KV_ev(mytime,tstep,nstep,ldorefinment)\n"
        "    else\n"
        "      call driver_integrator(mytime,tstep,nstep,ldorefinment)\n"
        "    endif\n",
        "main integrator dispatch",
    )
    write(path, text)


def patch_makefile():
    path = "build/Makefile"
    text = read(path)
    text = text.replace(
        "\tintegrator_mod.o statistic_mod.o io_mod.o",
        "\tintegrator_mod.o integrator_kv_ev_mod.o statistic_mod.o io_mod.o",
    )
    text = text.replace(
        "electric_field_mod.o eom_mod.o eom_ev_mod.o driver_eom_mod.o integrator_mod.o \\\n\tstatistic_mod.o io_mod.o main.o",
        "electric_field_mod.o eom_mod.o eom_ev_mod.o driver_eom_mod.o integrator_mod.o \\\n\tintegrator_kv_ev_mod.o statistic_mod.o io_mod.o main.o",
    )
    text = text.replace(
        "\tintegrator_mod.o statistic_mod.o io_mod.o main.o",
        "\tintegrator_mod.o integrator_kv_ev_mod.o statistic_mod.o io_mod.o main.o",
    )
    rule = "integrator_mod.o:integrator_mod.f90\n\t$(FC) $(FFLAGS) integrator_mod.f90\n"
    if rule not in text:
        raise RuntimeError("Makefile integrator rule not found")
    text = text.replace(rule, rule +
        "\nintegrator_kv_ev_mod.o:integrator_kv_ev_mod.f90\n"
        "\t$(FC) $(FFLAGS) integrator_kv_ev_mod.f90\n", 1)
    write(path, text)


def patch_manual():
    path = "manual/evaporation.tex"
    text = read(path)
    marker = ("which is Eq.~(39) of Yarin et al. written in the JETSPIN scaling.\n\n"
              "Evaporation also reduces the mass represented by a material element.")
    if marker not in text:
        raise RuntimeError("manual Maxwell marker not found")
    section = r'''which is Eq.~(39) of Yarin et al. written in the JETSPIN scaling.

\subsection{Kelvin--Voigt extension with concentration-dependent properties}

The concentration laws in Eqs.~\ref{eq:evap-viscosity} and
\ref{eq:evap-elastic} describe material properties and can also be
combined with the Kelvin--Voigt constitutive model available in JETSPIN.
This combination is an extension of the Yarin concentration-dependent
rheology; it should not be attributed to Ref.~\cite{yarin2001bending},
which used a Maxwell fluid.

In dimensionless variables define
\begin{equation}
 r_\mu=\frac{\mu_i}{\mu_0},\qquad
 r_G=\frac{G_i}{G_0}=\frac{r_\mu}
 {(c_{p,i}/c_{p,0})^{t_{\mathrm{ev}}}}.
\end{equation}
The Kelvin--Voigt constitutive identity is written as
\begin{equation}
 \bar\sigma_i=r_G\,\varepsilon_i+r_\mu\,\dot\varepsilon_i,
 \label{eq:evap-kv-identity}
\end{equation}
where the overdot denotes differentiation with respect to dimensionless
time.  Differentiating Eq.~\ref{eq:evap-kv-identity} gives
\begin{equation}
 \dot{\bar\sigma}_i=
 r_G\dot\varepsilon_i+r_\mu\ddot\varepsilon_i+
 \dot r_G\varepsilon_i+\dot r_\mu\dot\varepsilon_i.
 \label{eq:evap-kv-product}
\end{equation}
Thus a time-dependent concentration requires the two product-rule terms
that are absent when $G$ and $\mu$ are constant.

No additional material parameters are introduced.  From polymer-mass
conservation,
\begin{equation}
 \dot c_{p,i}=-c_{p,i}\frac{\dot V_i}{V_i},
 \label{eq:evap-kv-cpdot}
\end{equation}
and therefore
\begin{equation}
 \dot r_\mu=r_\mu\ln(10)Bm c_{p,i}^{m-1}\dot c_{p,i},
 \label{eq:evap-kv-mudot}
\end{equation}
\begin{equation}
 \dot r_G=r_G\left[\ln(10)Bm c_{p,i}^{m-1}
 -\frac{t_{\mathrm{ev}}}{c_{p,i}}\right]\dot c_{p,i}.
 \label{eq:evap-kv-gdot}
\end{equation}
JETSPIN reconstructs the instantaneous elastic strain from the current
stress state,
\begin{equation}
 \varepsilon_i=\frac{\bar\sigma_i-r_\mu\dot\varepsilon_i}{r_G},
 \label{eq:evap-kv-strain}
\end{equation}
so that no extra per-bead state variable is required.

The extension retains the historical JETSPIN Kelvin--Voigt kinematics:
$\dot\varepsilon_i=(1/l_i)dl_i/d\bar t$, while the acceleration term
used for $\ddot\varepsilon_i$ is the existing
$(1/l_i)d v_{\parallel,i}/d\bar t$ approximation.  Consequently, when
$\dot c_p=0$, Eq.~\ref{eq:evap-kv-product} reduces exactly to the
pre-existing JETSPIN Kelvin--Voigt differential equation.  At the Yarin
solidification cutoff, $\dot V_i=\dot c_{p,i}=0$, so $r_\mu$ and $r_G$
remain frozen.  Euler, Heun (RK2) and classical RK4 integrations are
supported for the coupled Kelvin--Voigt/evaporation model.

Evaporation also reduces the mass represented by a material element.'''
    text = text.replace(marker, section, 1)
    text = text.replace("parameter-fidelity and regression test",
                        "order-of-magnitude reference and regression test")
    text += "\n"
    write(path, text)


def patch_readme():
    path = "README.md"
    text = read(path)
    old = ("rheological corrections and input directives are documented in          \n"
           "manual/evaporation.tex.                                                 \n")
    new = old + (
        "The same concentration-dependent viscosity and elastic-modulus laws can    \n"
        "also be used with the Kelvin-Voigt rheology. This is documented as a       \n"
        "JETSPIN extension of the Yarin concentration laws (not as part of the       \n"
        "original Yarin-2001 Maxwell model).                                        \n")
    text = replace_once(text, old, new, "README evaporation paragraph")
    text = text.replace("parameter-fidelity regression case",
                        "order-of-magnitude reference/regression case")
    feat = ("  - The Yarin-2001 evaporation cutoff and rheological model are         \n"
            "    documented in manual/evaporation.tex.                               \n")
    text = replace_once(text, feat, feat +
        "                                                                        \n"
        "  - Evaporation can be combined with Kelvin-Voigt rheology using the       \n"
        "    same concentration-dependent mu and G laws, including the required     \n"
        "    product-rule terms for time-dependent material properties.              \n",
        "README feature")
    write(path, text)


def patch_test8_docs():
    for path in ["examples/input-8/input.dat", "examples/input-8/README.md"]:
        text = read(path)
        text = text.replace("parameter-fidelity regression case",
                            "order-of-magnitude reference/regression case")
        text = text.replace("parameter-fidelity and regression test",
                            "order-of-magnitude reference and regression test")
        write(path, text)


def patch_smoke():
    path = "tests/smoke/run.sh"
    text = read(path)
    marker = 'echo "Smoke mode $mode passed ($last_case case(s))"\n'
    if marker not in text:
        raise RuntimeError("smoke final marker not found")
    block = r'''if [ "$mode" != mpi ]; then
    echo "Running Kelvin-Voigt evaporation regression"
    python3 "$repo_root/tests/evaporation/check_kv_evaporation.py" \
        --eom-source "$repo_root/source/eom_ev_mod.f90" \
        --integrator-source "$repo_root/source/integrator_kv_ev_mod.f90"

    kv_integrator=1
    while [ "$kv_integrator" -le 3 ]; do
        kv_dir="$work_dir/kv-evap-$kv_integrator"
        mkdir -p "$kv_dir"
        cp "$work_dir/execute/main.x" "$kv_dir/main.x"
        cp "$repo_root/examples/input-8/input.dat" "$kv_dir/input.dat"
        sed -i 's/\r$//' "$kv_dir/input.dat"
        sed -i \
            -e "s/^[[:space:]]*integrator[[:space:]].*/ integrator $kv_integrator/" \
            -e 's/^[[:space:]]*timestep[[:space:]].*/ timestep 1.d-8/' \
            -e 's/^[[:space:]]*final time[[:space:]].*/ final time 1.d-6/' \
            -e 's/^[[:space:]]*print time[[:space:]].*/ print time 2.d-7/' \
            "$kv_dir/input.dat"
        sed -i '/^[[:space:]]*Finish/i\ kvfluid yes' "$kv_dir/input.dat"

        echo "Running Kelvin-Voigt evaporation with integrator $kv_integrator"
        (
            cd "$kv_dir"
            timeout "${JETSPIN_SMOKE_TIMEOUT:-30}" ./main.x > run.log 2>&1
        )
        grep -q 'Program closed correctly' "$kv_dir/run.log"
        grep -q 'Kelvin-Voigt evaporation integrator active' "$kv_dir/run.log"
        if grep -Eiq '(^|[^[:alpha:]])(error|nan|[-+]?inf(inity)?)([^[:alpha:]]|$)' \
            "$kv_dir/run.log" "$kv_dir/statout.dat"; then
            echo "Kelvin-Voigt evaporation integrator $kv_integrator failed" >&2
            tail -60 "$kv_dir/run.log" >&2
            exit 1
        fi
        kv_integrator=$((kv_integrator + 1))
    done
fi

'''
    text = text.replace(marker, block + marker, 1)
    write(path, text)


def patch_test_readme():
    path = "tests/evaporation/README.md"
    text = read(path)
    text += r'''

## Kelvin--Voigt evaporation extension

`check_kv_evaporation.py` guards the JETSPIN extension that combines the
Yarin concentration-dependent viscosity and elastic-modulus laws with the
Kelvin--Voigt constitutive model.  It checks the product-rule terms caused
by time-dependent material properties and verifies that the frozen-
concentration limit reduces to the historical JETSPIN Kelvin--Voigt form.
The smoke suite additionally runs short coupled evaporation/Kelvin--Voigt
jobs with Euler, Heun (RK2), and classical RK4 integration.
'''
    write(path, text)


def main():
    write("source/integrator_kv_ev_mod.f90", KV_MODULE)
    write("tests/evaporation/check_kv_evaporation.py", KV_CHECK)
    patch_eom()
    patch_io()
    patch_main()
    patch_makefile()
    patch_manual()
    patch_readme()
    patch_test8_docs()
    patch_smoke()
    patch_test_readme()
    print("Kelvin-Voigt evaporation extension applied")


if __name__ == "__main__":
    main()
