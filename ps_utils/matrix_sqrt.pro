function matrix_sqrt, matrix
  tol = 1e-8
  dblerr = 1e-10
  dim = size(matrix,/dimensions)
  if dim[0] ne dim[1] then begin
    message,'Matrix must be square'
  endif
  
  ; Calculate SVD of input matrix
  la_svd,matrix,svmatrix,wmatrix,vmatrix
   
  num = n_elements(svmatrix)
  
  ; verify that left and right singular vectors are equal
  svneg = dblarr(num)+1.0
  for ii = 0,num-1 do begin
    nerr = norm(wmatrix[ii,*]-vmatrix[ii,*])
    if nerr GT dblerr then svneg[ii] = -1D
  endfor
  negfind = where(svneg eq -1)
  if negfind[0] ne -1 then message,'Warning: Matrix has negative singular values and/or is non-symmetric'
  
  ; Scale Singular matrix with sqrt of singular values
  ;; dsvmatrix = dblarr(dim[0],dim[1])
  ;; dsvmatrix[indgen(dim[0]),indgen(dim[0])] = sqrt(svmatrix)
  dsvmatrix = diag_matrix(sqrt(svmatrix))
  sqmatrix = vmatrix##(dsvmatrix##la_invert(vmatrix))

   error = norm(sqmatrix##transpose(sqmatrix)-matrix)/norm(matrix)
  if (error ge tol) then message, 'Warning: ||LL^T-A||/||A|| > tol, tol: ' + string(tol) + ' error:' + string(error)
  return,sqmatrix
end