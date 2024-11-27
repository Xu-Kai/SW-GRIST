subroutine EXP(omnicopy)EXP(_)EXP(DTYPE)EXP(KIND)(dst, src, n)
  DTYPE(kind=KIND) :: dst(:), src(:)
  integer(kind=4) :: n
  call omnicopy_c(dst, src, n*KIND)
end subroutine EXP(omnicopy)EXP(_)EXP(DTYPE)EXP(KIND)
subroutine EXP(omnicopy)EXP(_s_)EXP(DTYPE)EXP(KIND)(dst, src)
  DTYPE(kind=KIND) :: dst(:), src(:)
  call omnicopy_c(dst, src, size(dst)*KIND)
end subroutine EXP(omnicopy)EXP(_s_)EXP(DTYPE)EXP(KIND)