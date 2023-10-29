     Program Sm
        INTEGER:: NI,NJ,MAXITER,I,J,key,NIO,NJO,N
        REAL(8):: Ls,hg,b,E,Br
        REAL(8):: Kappa,dx,dy,m,F,res,res1,dTheta1
        REAL(8),ALLOCATABLE,DIMENSION(:,:):: X,Y,P,Pn
        open(1,FILE='Input.txt')
        read(1,*) Ls
        read(1,*) hg
        read(1,*) b
        read(1,*) NI
        read(1,*) NJ
        read(1,*) Kappa
        read(1,*) m
        read(1,*) E
        read(1,*) MAXITER
        read(1,*) Br
        read(1,*) key
        close(1)

        allocate(X(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
        allocate(P(NI,NJ)) ! Pressure n-1
        allocate(Pn(NI,NJ))! Pressure n
!----------------- Coordinate of nodes ------------------------
        dx = (1.0d0*LS)/((NI-1))
        dy = 1.0d0/((NJ-1))
        DO I=1,NI
          DO J=1,NJ
            X(I,J) = (I-1)*dx
            Y(I,J) = (J-1)*dy
          END DO
        END DO
!---------- Coordinate of groove and holes ------------------------
        NIO = NINT(NI*((LS - b)/LS))
        NJO = NINT(NJ*(1.0d0 - b))
!----------------- Parameters ------------------------
        write(*,*)'W/L=',Ls,'h=',hg
        write(*,*)'NJ=',NJ,'NI=',NI,'NJO=',NJO,'NIO=',NIO
        write(*,*)'dx=',dx,'dy=',dy
!----------------- Initial fields -----------------------------
        P(1:NI,1:NJ) = 0.0d0
        Pn(1:NI,1:NJ) = 0.0d0
        call Border_conditions(NI,NJ,key,NIO,NJO,P)
!---------------- Solve equations ------------------------
        write(*,*) '-------Solve equations-------'
        Open(4,FILE='Res.plt')

        do N=1, MAXITER!метод простых итераций
            call Holes(NI,NJ,P,Pn,NIO,NJO,m,dx,dy,Kappa,key,dTheta1)!уравнение для отверстия
            call Groove(NI,NJ,P,Pn,NIO,NJO,Kappa,dx,dy,key,hg)!уравнение для канавок
            call Reynolds(NI,NJ,P,NIO,NJO,dx,dy,Pn)!уравнение Рейнольдса для расчетной области
            call Border_conditions(NI,NJ,key,NIO,NJO,Pn)!ГУ

            IF (N==1) THEN
                res1=MAXVAL(abs(Pn(:,:) - P(:,:)))
            ENDIF
            res = MAXVAL(abs(Pn(:,:) - P(:,:)))/res1
            CALL Force(NI,NJ,P,dx,dy,F)
            write(4,*) N, res, F
            !if (MOD(N,10000)==0) then
            !    write(*,*) N, res, P(1,NJO), P(NIO,1)
            !end if
            If (N/=MAXITER) then
                P(1:NI,1:NJ)=(1-Br)*P(1:NI,1:NJ)+Br*Pn(1:NI,1:NJ)
            end if
            if (res<E) then
                exit
            end if
        end do
        close(4)
        WRITE(*,*)'F=', F !несущая способность
        WRITE(*,*)'m=', m !коэф режима
 !----------------- Output data ------------------------------
        write(*,*) 'Output data'
        Open(2,FILE='Results.plt')
        Call Output_Fields(NI,NJ,X,Y,P,Pn)
        Close(2)
    END PROGRAM
!**********************************************************************************************
        SUBROUTINE Force(NI,NJ,P,dx,dy,F)
         INTEGER:: NI,NJ,I,J
         REAL(8):: dx,dy,F
         REAL(8),DIMENSION(NI,NJ):: P
         F = 0.0d0
         DO I=1,NI-1
             DO J=1,NJ-1
                 F=F+((P(I+1,J)+P(I,J)+P(I+1,J+1)+P(I,J+1))/4.0)*dx*dy
             END DO
         END DO
        END  SUBROUTINE
!**********************************************************************************************
        SUBROUTINE Output_Fields(NI,NJ,X,Y,P,Pn)
         INTEGER:: NI,NJ
         REAL(8),DIMENSION(NI,NJ):: X,Y,P,Pn
         Write(2,*) 'VARIABLES = "X", "Y", "P"'
         Write(2,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
         Write(2,'(100E25.16)') X(1:NI,1:NJ)
         Write(2,'(100E25.16)') Y(1:NI,1:NJ)
         Write(2,'(100E25.16)') P(1:NI,1:NJ)
        END  SUBROUTINE
!**********************************************************************************************
        SUBROUTINE Border_conditions(NI,NJ,key,NIO,NJO,P)
         INTEGER:: NI,NJ,I,J,key,NIO,NJO
         REAL(8),DIMENSION(NI,NJ):: P

         do I=1,NI
            P(I,NJ) = 0.0d0
         end do
         do J=1,NJ
            P(NI,J) = 0.0d0
         end do

         select case(key)
         case(1) !первый порядок точности гу
         do I=2,NIO-1
            P(I,1) = P(I,2)
         end do
         do I=NIO+1,NI-1
            P(I,1) = P(I,2)
         end do

         do J=1,NJO-1
            P(1,J) = P(2,J)
         end do
         do J=NJO+1,NJ-1
            P(1,J) = P(2,J)
         end do

         case(2) ! второй порядок точности гу
         do I=2,NIO-1
            P(I,1) = (4*P(I,2)-P(I,3))/3
         end do
         do I=NIO+1,NI-1
            P(I,1) = (4*P(I,2)-P(I,3))/3
         end do

         do J=1,NJO-1
            P(1,J) = (4*P(2,J)-P(3,J))/3
         end do
         do J=NJO+1,NJ-1
            P(1,J) = (4*P(2,J)-P(3,J))/3
         end do
         end select
        END  SUBROUTINE
!**********************************************************************************************
        subroutine Reynolds(NI,NJ,P,NIO,NJO,dx,dy,Pn)
         INTEGER:: NI,NJ,I,J,NIO,NJO
         REAL(8):: dx,dy
         REAL(8),DIMENSION(1:NI,1:NJ):: P,Pn

            do I=2,NIO-1
                do J=2,NJO-1
                    Pn(I,J)=(1/(2/(dx**2) + 2/(dy**2)))*(((P(I-1,J)+P(I+1,J))/(dx**2))+((P(I,J-1)+P(I,J+1))/(dy**2)))
                end do
            end do
            do I=NIO+1,NI-1
                do J=2,NJ-1
                    Pn(I,J)=(1/(2/(dx**2) + 2/(dy**2)))*(((P(I-1,J)+P(I+1,J))/(dx**2))+((P(I,J-1)+P(I,J+1))/(dy**2)))
                end do
            end do
            do I=2,NIO
                do J=NJO+1,NJ-1
                    Pn(I,J)=(1/(2/(dx**2) + 2/(dy**2)))*(((P(I-1,J)+P(I+1,J))/(dx**2))+((P(I,J-1)+P(I,J+1))/(dy**2)))
                end do
            end do
        end subroutine
 !**********************************************************************************************
        subroutine Holes(NI,NJ,P,Pn,NIO,NJO,m,dx,dy,Kappa,key,dTheta1)
         INTEGER:: NI,NJ,NIO,NJO,key
         REAL(8):: m,dx,dy,Theta1,Theta2,Kappa,dTheta1,dTheta2
         REAL(8),DIMENSION(NI,NJ):: P,Pn

         Theta1 = sqrt(1-P(1,NJO))!отверстие слева dx
         dTheta1 = abs(Theta1/(1-P(1,NJO)))
         Theta2 = sqrt(1-P(NIO,1))!отверстие справа dy
         dTheta2 = abs(Theta2/(1-P(NIO,1)))

         select case(key)
         case(1)
         Pn(1,NJO)=(1/(Kappa/dx + m*dTheta1))*(m*Theta1 + Kappa*P(2,NJO)/dx + m*P(1,NJO)*dTheta1)
         Pn(NIO,1)=(1/(Kappa/dy + m*dTheta2))*(m*Theta2 + Kappa*P(NIO,2)/dy + m*P(NIO,1)*dTheta2)

         case(2)
         Pn(1,NJO)=(m*Theta1 +m*dTheta1*P(1,NJO)+(0.5*Kappa/dx)*(4*P(2,NJO)-P(3,NJO)))*(1/(1.5*Kappa/dx + m*dTheta1))
         Pn(NIO,1)=(m*Theta2 +m*dTheta2*P(NIO,1)+(0.5*Kappa/dy)*(4*P(NIO,2)-P(NIO,3)))*(1/(1.5*Kappa/dy + m*dTheta2))
         end select
        end subroutine
 !**********************************************************************************************
        subroutine Groove(NI,NJ,P,Pn,NIO,NJO,Kappa,dx,dy,key,h)
         INTEGER:: NI,NJ,I,J,NIO,NJO,key
         REAL(8):: Kappa,dx,dy,h
         REAL(8),DIMENSION(NI,NJ):: P,Pn

         selectcase(key)
         case(1)
         do I=2,NIO-1!верхняя канавка
            Pn(I,NJO)=(-((h**3)/dy)*(P(I,NJO+1)+P(I,NJO-1))-(Kappa/(dx**2))*(P(I+1,NJO)+P(I-1,NJO)))/(2*(-((h**3)/dy)-(Kappa/(dx**2))))
         end do
         do J=2,NJO-1 !нижняя канавка
            Pn(NIO,J)=(-((h**3)/dx)*(P(NIO+1,J)+P(NIO-1,J))-(Kappa/(dy**2))*(P(NIO,J-1)+P(NIO,J+1)))/(2*(-((h**3)/dx)-(Kappa/(dy**2))))
         end do
         Pn(NIO,NJO)=(P(NIO-1,NJO)/dx +P(NIO,NJO-1)/dy)/(1/dx +1/dy)

         case(2)
         do I=2,NIO-1!верхняя канавка
            Pn(I,NJO)=(-(P(I-1,NJO)+P(I+1,NJO))+(h*h*h*dx*dx/(2*Kappa*dy))*(-4*P(I,NJO-1)+P(I,NJO-2)-4*P(I,NJO+1)+P(I,NJO+2)))/(-2-3*h*h*h*dx*dx/(Kappa*dy))
         end do
         do J=2,NJO-1 !нижняя канавка
            Pn(NIO,J)=(-(P(NIO,J-1)+P(NIO,J+1))+(h*h*h*dy*dy/(2*Kappa*dx))*(-4*P(NIO-1,J)+P(NIO-2,J)-4*P(NIO+1,J)+P(NIO+2,J)))/(-2-3*h*h*h*dy*dy/(Kappa*dx))
         end do
         Pn(NIO,NJO)=((4*P(NIO-1,NJO)-P(NIO-2,NJO))/dx +(4*P(NIO,NJO-1)-P(NIO,NJO-2))/dy)/(3/dx +3/dy)
         end select
        end subroutine


