 module parse_module
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which parse the input file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
 private
 
 public :: getline,intstr,dblstr,strip,lowcase,copystring,getword
 public :: findstring,checkempty,findwords
 
 contains

 subroutine getline(safe,ifile,lenstring,string)
 
!***********************************************************************
!     
!     JETSPIN subroutine for reading a text line and store it as
!     character string
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
      
  logical, intent(out) :: safe
  integer, intent(in) :: ifile
  integer, intent(in) :: lenstring
  character(len=lenstring), intent(out) :: string
  
  character(len=6) :: string6  
      
  safe=.true.
  
  write(string6,'(a,i3,a)')'(a',lenstring,')'
  read(ifile,string6,end=100)string
  
  return
  
 100  safe=.false.
 
  return
       
 end subroutine getline

 function intstr(string,lenstring,laststring)
 
!***********************************************************************
!     
!     JETSPIN subroutine for extracting integers from a character 
!     string
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(inout) :: string
  integer, intent(in) :: lenstring
  integer, intent(out) :: laststring
  
  integer :: intstr
  
  integer :: j,isn
  character*1, parameter, dimension(0:9) :: & 
   n=(/'0','1','2','3','4','5','6','7','8','9'/)
  logical :: flag,lcount,final
  character*1 :: ksn
  character*1, dimension(lenstring) :: word
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo
  
  isn=1
  laststring=0
  ksn='+'
  intstr=0
  flag=.false.
  final=.false.
  lcount=.false.
  
  
  do while(laststring<lenstring.and.(.not.final))
    
    laststring=laststring+1
    flag=.false.
    
    do j=0,9
      
      if(n(j)==word(laststring))then
        
        intstr=10*intstr+j
        lcount=.true.
        flag=.true.
        
      endif
    
    enddo
    
    if(lcount.and.(.not.flag))final=.true.
    if(flag .and. ksn=='-')isn=-1
    ksn=word(laststring)
    
  enddo

  intstr=isn*intstr

  do j=laststring,lenstring
    word(j-laststring+1)=word(j)
  enddo
  do j=lenstring-laststring+2,lenstring
    word(j)=' '
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
  end function intstr

  function dblstr(string,lenstring,laststring)
  
!***********************************************************************
!     
!     JETSPIN subroutine for extracting double precisions from a  
!     character string
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(inout) :: string
  integer, intent(in) :: lenstring
  integer, intent(out) :: laststring
  
  double precision :: dblstr
  
  logical :: flag,ldot,start,final
  integer :: iexp,idum,i,j,fail
  double precision :: sn,ten,one

  character*1, parameter, dimension(0:9) :: & 
   n=(/'0','1','2','3','4','5','6','7','8','9'/)
  character*1, parameter :: dot='.'
  character*1, parameter :: d='d'
  character*1, parameter :: e='e'
  
  character*1 :: ksn
  character*1, dimension(lenstring) :: word
  character(len=lenstring) :: work
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo

  laststring=0
  sn=1.d0
  ksn='+'
  ten=10.d0
  one=1.d0
 
  dblstr=0.d0
  iexp=0
  idum=0
  start=.false.
  ldot=.false.
  final=.false.
  
  do while(laststring<lenstring .and. (.not.final))
    
    laststring=laststring+1
    flag=.false.
    
    do j=0,9
      
      if(n(j)==word(laststring))then
        
        dblstr=ten*dblstr+one*dble(j)
        flag=.true.
        start=.true.
            
      endif
          
    enddo
        
    if(dot==word(laststring))then
          
      flag=.true.
      ten=1.d0
      ldot=.true.
      start=.true.
          
    endif

    if(flag .and. ksn=='-') sn=-1.d0
    if(ldot)one=one/10.d0
    ksn=word(laststring)
    if(ksn=="D")ksn="d"
    if(ksn=="E")ksn="e"
    
    if(start)then
      if(d==ksn .or. e==ksn)then
        do i=1,lenstring-laststring
          work(i:i)=word(i+laststring)
        enddo
        iexp=intstr(work,lenstring-laststring,idum)
        final=.true.
      endif
      if(.not.flag)final=.true.        
    endif
  enddo
  
  dblstr=sn*dblstr*(10.d0**iexp)
  laststring=laststring+idum
  
  do j=laststring,lenstring
    word(j-laststring+1)=word(j)
  enddo
  do j=lenstring-laststring+2,lenstring
    word(j)=' '
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
 end function dblstr

 subroutine strip(string,lenstring)

  implicit none
  
  character(len=*) :: string
  integer, intent(in) :: lenstring
  
  integer :: i,j
  character*1, dimension(lenstring) :: word
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo 
  
  do i=1,lenstring
    
    if(word(1)==' ')then
      
      do j=1,lenstring-1
        
        word(j)=word(j+1)
        
      enddo
      
      word(lenstring)=' '
      
    endif
    
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
 end subroutine strip

 subroutine lowcase(string,lenstring)
 
!***********************************************************************
!     
!     JETSPIN subroutine to lowercase a character string
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(inout) :: string
  integer, intent(in) :: lenstring
  
  character*1, dimension(lenstring) :: word
  character*1 :: letter
  
  integer :: i,j
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo 
  
  do i=1,lenstring
    
    letter=word(i)
    
    if(letter=='A')then
      letter='a'
    else if(letter=='B')then
      letter='b'
    else if(letter=='C')then
      letter='c'
    else if(letter=='D')then
      letter='d'
    else if(letter=='E')then
      letter='e'
    else if(letter=='F')then
      letter='f'
    else if(letter=='G')then
      letter='g'
    else if(letter=='H')then
      letter='h'
    else if(letter=='I')then
      letter='i'
    else if(letter=='J')then
      letter='j'
    else if(letter=='K')then
      letter='k'
    else if(letter=='L')then
      letter='l'
    else if(letter=='M')then
      letter='m'
    else if(letter=='N')then
      letter='n'
    else if(letter=='O')then
      letter='o'
    else if(letter=='P')then
      letter='p'
    else if(letter=='Q')then
      letter='q'
    else if(letter=='R')then
      letter='r'
    else if(letter=='S')then
      letter='s'
    else if(letter=='T')then
      letter='t'
    else if(letter=='U')then
      letter='u'
    else if(letter=='V')then
      letter='v'
    else if(letter=='W')then
      letter='w'
    else if(letter=='X')then
      letter='x'
    else if(letter=='Y')then
      letter='y'
    else if(letter=='Z')then
      letter='z'
    endif
    
    word(i)=letter
    
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
 end subroutine lowcase
 
 subroutine copystring(oldstring,newstring,lenstring)
 
!***********************************************************************
!     
!     JETSPIN subroutine to copy one character string into another
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(in) :: oldstring
  character(len=*), intent(out) :: newstring
  integer, intent(in) :: lenstring
  
  integer :: i
  
  do i=1,lenstring
    newstring(i:i)=oldstring(i:i)
  enddo
  
  return
  
 end subroutine copystring
  
 function findstring(seek,string,here,lenstring)
 
!***********************************************************************
!     
!     JETSPIN subroutine to find an explicit string in an input record
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(in) :: seek
  character(len=*), intent(inout) :: string
  integer, intent(out) :: here
  integer, intent(in) :: lenstring
  
  logical :: findstring
  character*1, dimension(lenstring) :: word
  integer :: i,j,nseek,m
  logical :: findspace
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo 

  m=lenstring
  nseek=len(seek)
  findstring=.false.
  
  here=0
  findspace=.true.
  do while(here<m-nseek .and. (.not.findstring))
    
    findstring=.true.
    
    do i=1,nseek
      if(seek(i:i)/=word(here+i))findstring=.false.
    enddo
    findstring=(findstring.and.findspace)

    here=here+1
    
    if(word(here)==' ')then
      findspace=.true.
    else
      findspace=.false.
    endif

  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo

  return
  
 end function findstring

 subroutine striptext(string,lenstring,nwords)
 
!***********************************************************************
!     
!     JETSPIN subroutine to strip leading text from a data record
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(inout) :: string
  integer, intent(in) :: lenstring,nwords
  
  logical :: final
  integer :: i,j,k
  character*1, dimension(lenstring) :: word

  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo
  
  do k=1,nwords
  
    i=0
    final=.false.
    
    do while((.not.final) .and. i<lenstring)
      
      i=i+1
      
      if(word(1)==' ')then
        
        final=.true.
        
      else
      
        do j=1,lenstring-1
          
          word(j)=word(j+1)
          
        enddo
        
        word(lenstring)=' '
        
      endif
      
    enddo
    
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
  end subroutine striptext

 subroutine getword(stringout,string,len1,len2)
 
!***********************************************************************
!     
!     JETSPIN subroutine to fetch a character word from a string
!     while ignoring leading blanks
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: len1,len2
  character(len=len1), intent(out) :: stringout
  character(len=*), intent(in) :: string
  
  
  character*1, dimension(len1) :: wrdseq
  character*1, dimension(len2) :: word
  
  logical :: final
  integer :: i,j,k
  
  do j=1,len2
    word(j)=string(j:j)
  enddo
  
  do i=1,len1
    wrdseq(i)=' '
  enddo
  
  i=0
  k=0
  final=.false.
  
  do while((.not.final) .and. i<len2)
    
    i=i+1
    
    if(word(1)==' ')then
      
      if(k>0)final=.true.
      
    else
      
      k=k+1
      wrdseq(k)=word(1)
      if(k==len1)final=.true.
      
    endif
    
    do j=1,len2-1
      
      word(j)=word(j+1)
      
    enddo
    
    word(len2)=' '
    
  enddo
  
  do j=1,len1
    stringout(j:j)=wrdseq(j)
  enddo

  return
  
 end subroutine getword
 
 subroutine findwords(nwords,outwords,string,lenstring)
 
!***********************************************************************
!     
!     JETSPIN subroutine to count the number of words contained in a
!     character string and return them as elements of a character array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: lenstring
  character(len=*), intent(in) :: string
  integer, intent(out) :: nwords
  character(len=lenstring) ,allocatable, intent(inout) :: outwords(:)
  
  character(len=lenstring) :: temps,temps2
  character*1, allocatable, dimension(:,:) :: tempswords,foundwords
  integer :: i,j,k
  logical :: lredo
  character(len=3) :: temps3
  character(len=7) :: temps7
  
  temps(1:lenstring)=string(1:lenstring)
  
  i=0
  do while(.not.checkempty(temps))
    i=i+1
    if(allocated(foundwords))then
      allocate(tempswords(1:lenstring,1:i-1))
      tempswords(1:lenstring,1:i-1)=foundwords(1:lenstring,1:i-1)
      deallocate(foundwords)
      allocate(foundwords(1:lenstring,1:i))
      foundwords(1:lenstring,1:i-1)=tempswords(1:lenstring,1:i-1)
      deallocate(tempswords)
    else
      allocate(foundwords(1:lenstring,1:i))
    endif
    temps=adjustl(temps)
    j=0
    lredo=.true.
    do while(lredo)
      j=j+1
      if(temps(j:j).ne.' ')then
        foundwords(j,i)=temps(j:j)
      else
       do k=j,lenstring
         foundwords(k,i)=' '
       enddo
       do k=1,lenstring-j
         temps2(k:k)=temps(j+k-1:j+k-1)
       enddo
       do k=lenstring-j+1,lenstring
         temps2(k:k)=' '
       enddo
       temps(:)=temps2(:)
       lredo=.false.
      endif
    enddo
  enddo
  
  nwords=i
  
  if(allocated(outwords))then
    deallocate(outwords)
  endif
  allocate(outwords(1:nwords))
  
  write(temps3,'(i3)')lenstring
  write(temps7,'(a1,a3,a3)')'(',adjustr(temps3),'a1)'
  
  do i=1,nwords
    write(temps,temps7)(foundwords(j,i),j=1,lenstring)
    outwords(i)=temps
  enddo
  
  return
  
 end subroutine findwords
 
 function checkempty(string)
 
!***********************************************************************
!     
!     JETSPIN function for checking if a character string is empty
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(in) :: string
  
  logical :: checkempty
  
  integer :: i
  
  i=len(trim(string))
  
  checkempty=(i==0)

  return
  
 end function checkempty
 
 end module parse_module

