C FILE: FIB1.F
      SUBROUTINE FIB(A,N)
C
C     CALCULATE FIRST N FIBONACCI NUMBERS
C
      INTEGER N
      REAL*8 A(N)
      DO I=1,N
         IF (I.EQ.1) THEN
            A(I) = 0.0D0
         ELSEIF (I.EQ.2) THEN
            A(I) = 1.0D0
         ELSE 
            A(I) = A(I-1) + A(I-2)
         ENDIF
      ENDDO
      END


      SUBROUTINE make2(g1,gmv,coeff,size_u,size_a,size_i,g2,g2mv)
      IMPLICIT NONE      
      INTEGER size_a, size_u, size_i,u,v,w,x,m,n,o
      REAL*8:: mat_ele=0.0d0,ele=0.0d0
      REAL*8,DIMENSION(:)::coeff
      REAL*8, DIMENSION(:,:) ::g1
!      REAL*8, DIMENSION(:,:,:,:) ::g2

      REAL*8 :: g2(size_u,size_u,size_u,size_u)
      REAL*8 :: g2mv(size_u,size_u,size_u,size_u,size_a,size_a)
      REAL*8, DIMENSION(:,:,:,:) ::gmv
!      REAL*8, DIMENSION(:,:,:,:,:,:) ::g2mv

Cf2py intent(out) g2
Cf2py intent(out) g2mv
      write(*,*) "making gamma2"

!      write (*,*) coeff(1)
      !write (*,*) g1(1,2)
      
      do u =1,size_u
       do v =1,size_u
        do w =1,size_u
         do x =1,size_u
          do m =1,size_a
           do n =1,size_a
            do o =1,size_a
             mat_ele=mat_ele+coeff(m)*coeff(n)*gmv(u,x,m,o)*gmv(v,w,o,n)
            end do
            if (v == x) then
             mat_ele=mat_ele-coeff(m)*coeff(n)*gmv(u,w,m,n)
            end if
            !if (mat_ele>0.000001d0) then
             g2mv(u,v,w,x,m,n) = mat_ele
             ele=ele+mat_ele
             mat_ele = 0.0d0
            !else
            ! g2mv(u,v,w,x,m,n)=0.d0
            ! mat_ele=0.d0
            !end if
           end do
          end do
          !write(*,*) u,v,w,x,ele
          !if (ele>0.000001d0) then
           g2(u,v,w,x) = ele
          !else
          ! g2(u,v,w,x)=0.d0
          !end if
          ele = 0.0d0
         end do
        end do
       end do
      end do
      END

      SUBROUTINE lambda2(g1,g2, coeff, size_u, size_a, size_i, l2)

      INTEGER :: u,v,w,x,size_u,size_i,size_a
      REAL*8, Dimension(:,:) :: g1
      REAL*8, Dimension(:) :: coeff
      REAL*8, Dimension(:,:,:,:) :: g2
      REAL*8:: ele=0.d0
      REAL*8 :: l2(size_u,size_u,size_u,size_u)
Cf2py intent(out) l2
      do u =1,size_u
       do v =1,size_u
        do w =1,size_u
         do x =1,size_u
          ele = 0.d0
          write(*,*) u,v,w,x,g2(u,v,w,x)
          ele = g2(u,v,w,x)-g1(u,x)*g1(v,w) 

          !write(*,*) "ele1 ",ele
          ele =ele + 1.0d0/2.0d0*g1(u,w)*g1(v,x)
          !write(*,*) "ele1+ele2 ",ele
          !if (ele>0.000001d0) then
           l2(u,v,w,x)=ele
          !else
           !l2(u,v,w,x)=0.d0
          !end if
         end do
        end do
       end do
      end do


      write(*,*) "check for 1: ", g2(1,1,1,2)-g1(1,2)*g1(1,1),
     +   g2(1,1,1,2),g1(1,2),g1(1,1) 
      write(*,*) "check for 1: ", 1.d0/2.d0*g1(1,2)*g1(1,1) 
      write(*,*) "check for 2: ", g2(2,1,1,1)-g1(2,1)*g1(1,1),
     +   g2(2,1,1,1),g1(1,2),g1(1,1) 
      write(*,*) "check for 2: ", 1.d0/2.d0*g1(2,1)*g1(1,1) 
      END
      SUBROUTINE make3(gmv,g2,g2mv,coeff,size_u,size_a,size_i,g3)
      IMPLICIT NONE      
      INTEGER size_a, size_u, size_i,u,v,w,x,y,z,m,n,o
      REAL*8:: mat_ele=0.0d0
      REAL*8,DIMENSION(:)::coeff
      REAL*8, DIMENSION(:,:,:,:) ::gmv
      REAL*8, DIMENSION(:,:,:,:) ::g2

      !REAL*8 :: g2(size_u,size_u,size_u,size_u)

      !REAL*8 :: g2mv(size_u,size_u,size_u,size_u,size_a,size_a)
      REAL*8 :: g3(size_u,size_u,size_u,size_u,size_u,size_u)
      REAL*8, DIMENSION(:,:,:,:,:,:) ::g2mv
!      REAL*8, DIMENSION(:,:,:,:,:,:) ::g2mv

Cf2py intent(out) g3
      write(*,*) "Gamma 3 values"

!      write (*,*) coeff(1)
      !write (*,*) g1(1,2)
      
      do u =1,size_u
       do v =1,size_u
        do w =1,size_u
         do x =1,size_u
          do y =1,size_u
           do z =1,size_u
            do m =1,size_a
             do n =1,size_a
              do o =1,size_a
              
c               mat_ele=mat_ele+coeff(m)*coeff(n)*g2mv(u,v,x,y,m,o)
c     +           *gmv(w,z,o,n)
c              end do
c               if (w == y) then
c                mat_ele=mat_ele-coeff(m)*coeff(n)*g2mv(u,v,x,z,m,n)
c               end if
c               if (w == x) then
c                mat_ele=mat_ele+coeff(m)*coeff(n)*g2mv(u,v,y,z,m,n)
c		'''
               mat_ele=mat_ele+coeff(m)*coeff(n)*g2mv(u,v,z,y,m,o)
     +           *gmv(w,x,o,n)
              end do
               if (w == y) then
                mat_ele=mat_ele-coeff(m)*coeff(n)*g2mv(u,v,z,x,m,n)
               end if
               if (w == z) then
                mat_ele=mat_ele+coeff(m)*coeff(n)*g2mv(u,v,x,y,m,n)
               end if
             end do
            end do
            !write(*,*) u,v,w,x,mat_ele
            g3(u,v,w,x,y,z) = mat_ele
            mat_ele = 0.0d0
           end do
          end do
         end do
        end do
       end do
      end do
      END

      SUBROUTINE lambda3(g1,g3,l2,coeff, size_u, size_a, size_i, l3)

      INTEGER :: u,v,w,x,y,z,size_u,size_i,size_a
      REAL*8, Dimension(:,:) :: g1
      REAL*8, Dimension(:) :: coeff
      REAL*8 :: ele = 0.d0
      REAL*8, Dimension(:,:,:,:,:,:) :: g3
      REAL*8, Dimension(:,:,:,:) :: l2

      REAL*8 :: l3(size_u,size_u,size_u,size_u,size_u,size_u)
Cf2py intent(out) l3
      do u =1,size_u
       do v =1,size_u
        do w =1,size_u
         do x =1,size_u
          do y =1,size_u
           do z =1,size_u
            ele = 0.d0
            ele=g3(u,v,w,x,y,z) - g1(u,x)*l2(v,w,y,z) - 
     +          g1(v,y)*l2(u,w,x,z) - g1(w,z)*l2(u,v,x,y)
            ele = ele+1.d0/2.d0*g1(u,y)*l2(v,w,x,z)+
     +     1.d0/2.d0*g1(u,z)*l2(v,w,x,y)+1.d0/2.d0*g1(v,x)*l2(u,w,y,z)+
     +     1.d0/2.d0*g1(v,z)*l2(u,w,x,y)+1.d0/2.d0*g1(w,x)*l2(u,v,y,z)+
     +     1.d0/2.d0*g1(w,y)*l2(u,v,x,z)

            ele=ele-g1(u,x)*g1(v,y)*g1(w,z)

            ele=ele+1.d0/2.d0*g1(u,x)*g1(v,z)*g1(w,y)+
     +     1.d0/2.d0*g1(u,z)*g1(v,y)*g1(w,x)+
     +     1.d0/2.d0*g1(u,y)*g1(v,x)*g1(w,z)

            ele=ele-1.d0/4.d0*g1(u,y)*g1(v,z)*g1(w,x)-
     +     1.d0/4.d0*g1(u,z)*g1(v,x)*g1(w,y)
            l3(u,v,w,x,y,z)=ele
           end do
          end do
         end do
        end do
       end do
      end do
      write(*,*) "l3 exits"
      END
      SUBROUTINE vector_comp(l3,g3,size_u,val)

      INTEGER :: u,v,w,x,y,z,size_u
      REAL*8 :: ele = 0.d0, val
      REAL*8, Dimension(:,:,:,:,:,:) :: l3
Cf2py intent(out) val
      do u =1,size_u
       do v =1,size_u
        do w =1,size_u
         do x =1,size_u
          do y =1,size_u
           do z =1,size_u
            !write(*,*) g3(u,v,w,z,y,x), g3(x,y,z,w,v,u)
            ele=ele+l3(u,v,w,z,y,x)*l3(x,y,z,w,v,u)
           end do
          end do
         end do
        end do
       end do
      end do
      val=sqrt(ele)
      write(*,*) "comparing value = ", val
      END
C END FILE FIB1.F
