! Convert an integer to a string
      subroutine int2str(inint,nchar,str)
      character*80 str

      itmp=abs(inint)
      str=''
      do i=1,nchar
         j=mod(itmp,10)
         str=achar(48+j)//str
         itmp=(itmp-j)/10
      end do

      if(inint.lt.0) str='-'//str

      return
      end
