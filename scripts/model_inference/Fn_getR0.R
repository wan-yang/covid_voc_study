
Fn_getR0_SEIR=function(PARMS){
  with(as.list(PARMS), {
    Ro = beta.mean * Tir.mean
    
    Ro
  })
}


Fn_getR0_SEIRageGr=function(PARMS){
  with(as.list(PARMS), {
    
    Ro.matrix = matrix(0,num_gr,num_gr)
    for(i in 1:num_gr){
      for(j in 1:num_gr){
        # r_ij: same famula for without age structure
        r_ij = BETA.I[i,j] * Tir.mean
        Ro.matrix[i,j]=Npops[i]/Npops[j]*r_ij
      }
    }
    Ro=eigen(Ro.matrix)$values[1]
    Re(Ro)
  })
}


Fn_getRe_SEIRageGr=function(PARMS){
  with(as.list(PARMS), {
    
    Ro.matrix = matrix(0,num_gr,num_gr)
    for(i in 1:num_gr){
      for(j in 1:num_gr){
        # r_ij: same famula for without age structure
        r_ij = BETA.I[i,j] * Tir.mean[j]  # is Tir.mean[j] or Tir.mean[i]?? - should be the infector, so Tir.mean[j]
        Ro.matrix[i,j]=S[i]/Npops[j]*r_ij  # replace Npops[i] with S[i]
      }
    }
    Ro=eigen(Ro.matrix)$values[1]
    Re(Ro)
  })
}
