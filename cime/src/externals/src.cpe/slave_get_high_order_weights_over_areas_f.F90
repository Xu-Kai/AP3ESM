    subroutine slave_cal(num_area,num_seg_max,irecons,ngpc,x,dx,gspts,weights,gsweights,num_seg,iarea)
    use shr_kind_mod,   only: r8=>shr_kind_r8
    implicit none
    
    integer ,intent(in):: num_area, num_seg_max, irecons, ngpc,iarea
    integer   , dimension(num_area), intent(in)   :: num_seg
    
    real(kind=r8), dimension(2,num_seg_max,num_area), intent(in):: x, dx
    real(kind=r8), dimension(ngpc), intent(in):: gsweights, gspts
    
    real(kind=r8), dimension(irecons,num_area), intent(out)  :: weights

    real(kind=r8), dimension(ngpc,num_seg_max               ) :: xq,yq
    real(kind=r8), dimension(ngpc,num_seg_max,irecons) :: F
    real(kind=r8), dimension(irecons) :: weights_area
    real(kind=r8), dimension(ngpc,num_seg_max) :: xq2, yrh, rho, tmp
    real(kind=r8) , dimension(ngpc,num_seg_max) :: xq2ir, xq2i, rhoi
    integer :: iseg,i,j,k
    

    do iseg=1,num_seg(iarea)
    !!write(*,*) 'iseg', iseg
        xq(:,iseg) = x(1,iseg,iarea)+dx(1,iseg,iarea)*gspts(:)
        yq(:,iseg) = x(2,iseg,iarea)+dx(2,iseg,iarea)*gspts(:)
    end do
    do j=1,num_seg(iarea)
        do i=1,ngpc
          xq2(i,j)   =  xq(i,j)*xq(i,j)
          xq2i(i,j)  =  1.0_r8/(1.0_r8+xq2(i,j))
          xq2ir(i,j) =  SQRT(xq2i(i,j))
          rho(i,j)   =  SQRT(1.0_r8+xq2(i,j)+yq(i,j)*yq(i,j))
          rhoi(i,j)  =  1.0_r8/rho(i,j)
          yrh(i,j)   =  yq(i,j)*rhoi(i,j)
          tmp(i,j)   =  yq(i,j)*xq2ir(i,j)
          F(i,j,1)   =  yrh(i,j)*xq2i(i,j)
          F(i,j,2)   =  xq(i,j)*yrh(i,j)*xq2i(i,j)
          F(i,j,3)   = -1.0_r8*rhoi(i,j)
          F(i,j,4)   =  xq2(i,j)*yrh(i,j)*xq2i(i,j)
          F(i,j,6)   = -xq(i,j)*rhoi(i,j)
        enddo

        do i=1,ngpc
            F(i,j,5)   = -yq(i,j)*rhoi(i,j)+log(tmp(i,j)+rho(i,j)*xq2ir(i,j))
        end do
    enddo

    weights_area = 0.0_r8

    do k=1,irecons
        do iseg=1,num_seg(iarea)
            weights_area(k) = weights_area(k) + sum(gsweights(:)*F(:,iseg,k))*0.5_r8*dx(1,iseg,iarea)
        end do
    end do
    weights(1:irecons,iarea) = weights_area(1:irecons)
    
    end subroutine slave_cal