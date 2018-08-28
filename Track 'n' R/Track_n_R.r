################################################################################
#                                                                              #
#                                  Track 'n' R                                 #
#                                                                              #
#                   Core functions linking PointTracker to R                   #
#                   for analysis of cell division and growth                   #
#                                                                              #
#                             Florent Pantin, 2013                             #
#                                                                              #
################################################################################



## Libraries ##

options(repos = "http://cran.ma.imperial.ac.uk/") # UK/London
for (pkg in c("spatial", "sp", "splancs", "rgdal", "RSEIS", "rgeos", "grDevices", "rJava"))
  {
  if (!pkg %in% installed.packages()[, "Package"]) { install.packages(pkg) }
  #update.packages(pkg)
  library(pkg, character.only = T)
  }



## Functions ##


################################################################################
#                                     Basic                                    #
################################################################################

# This function defines an operator allowing to retrieve more than one value
# from a function output without having to use 'unlist'.
# See http://code.google.com/p/miscell/source/browse/rvalues/rvalues.r
':=' = function(lhs, rhs)
  {
  frame = parent.frame()
  lhs = as.list(substitute(lhs))
  if (length(lhs) > 1)
  lhs = lhs[-1]
  if (length(lhs) == 1) {
  do.call(`=`, list(lhs[[1]], rhs), envir=frame)
  return(invisible(NULL)) }
  if (is.function(rhs) || is(rhs, 'formula'))
  rhs = list(rhs)
  if (length(lhs) > length(rhs))
  rhs = c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
  do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL))
  }



################################################################################
#          R translation of Pierre Barbier de Reuille's growth_algo.py         #
################################################################################


linearToExponential_growth <- function (value, DT)
  {
  # Correct a (list of) value(s) for growth estimated using linear-growth hypothesis
  # into value(s) corresponding to an exponential growth hypothesis.
  #
  # :Note: You should convert growth rate and vorticity but not the direction.
  # i.e. Do not convert directly the full tensor as it would not work.

  return (log(value*DT+1)/DT)
  }


exponentialToLinear_growth <- function (value, DT)
  {
  return ((exp(value*DT)-1)/DT)
  }


linearToExponential_strain <- function (tr, DT)
  {
  eigentr <- eigen(tr, symmetric = T)
  w <- eigentr$values
  v <- eigentr$vectors
  w <- linearToExponential_growth(w, DT)
  return (tr = v %*% diag(w) %*% t(v))
  }


exponentialToLinear_strain <- function (tr, DT)
  {
  eigentr <- eigen(tr, symmetric = T)
  w <- eigentr$values
  v <- eigentr$vectors
  w <- exponentialToLinear_growth(w, DT)
  return (tr = v %*% diag(w) %*% t(v))
  }


linearToExponential_tensor <- function (tr, DT)
  {
  s <- (tr+t(tr))/2
  v <- (tr-t(tr))/2
  return (nt = v + linearToExponential_strain(s, DT))
  }


exponentialToLinear_tensor <- function (tr, DT)
  {
  s <- (tr+t(tr))/2
  v <- (tr-t(tr))/2
  return (nt = v + exponentialToLinear_strain(s, DT))
  }


fitmat <- function (P, Q)
  {
  # Compute the best transformation to map P onto Q without translation.
  #
  # The transformation is best according to the least square residuals criteria and correspond to the matrix M in:
  #
  #      P-p0 = M(Q-q0) + T
  #
  # where p0 and q0 are the barycentre of the list of points P and Q.
  #
  # Parameters:
  #      P : array(M,2)
  #          List of points p, one point per line, first column for x, second for y
  #      Q : array(M,2)
  #          List of points q, one point per line, first column for x, second for y

  p0 <- matrix(colMeans(P), nrow = nrow(P), ncol = ncol(P), byrow = T)
  q0 <- matrix(colMeans(Q), nrow = nrow(Q), ncol = ncol(Q), byrow = T)
  P <- P - p0
  Q <- Q - q0
  A <- t(P) %*% P
  if (1/rcond(A, "1") > 1e15) { return ("None") }
  V <- t(P) %*% Q
  M <- (solve(A) %*% V)
  return (t.M = t(M))
  }


polarDecomposition <- function (M, at_start = T)
  {
  # Decompose the 2D/3D transformation matrix M (i.e. 2x2 or 3x3 matrix) into a rotation and a symmetric matrix.
  #
  # The symmetric matrix is such that only the smallest eigenvalue may be negative.
  #
  # :returns: (S,R), with S the scaling and R the rotation.
  # :returntype: (matrix of float,matrix of float)
  #
  #  If at_start is True, then the scaling is performed first (i.e. M = R*S).

  SVD <- svd(M)
  w <- SVD$u
  v <- SVD$v
  d <- SVD$d
  R <- w %*% v
  if (det(R) < 0) # Is the unit matrix not a rotation (i.e. direct)
    {
    d[length(d)] <- -d[length(d)]
    v[, ncol(v)] <- -v[, ncol(v)]
    R <- w %*% v
    }
  d <- diag(d)
  if (at_start) { S <- Conj(t(v)) %*% d %*% v }
  else { S <- w %*% d %*% Conj(t(w)) }
  return (list(S = S, R = R))
  }


rotationToVorticity <- function (R, DT)
  {
  if (ncol(R) == 1 & nrow(R) == 1)
    {
    return (0) # there is no rotation ...
    }
  else if (ncol(R) == 2 & nrow(R) == 2)
    {
    A <- atan2(R[2,1], R[1,1])
    R <- matrix(c(0, -A, A, 0), ncol = 2, byrow = T)
    return (list(RoverDT = R/DT, A = A))
    }
  else if (ncol(R) == 3 & nrow(R) == 3)   ### Warning: not tested ###
    {
    eigenR <- eigen(R, symmetric = F)
    ev <- eigenR$values
    ec <- eigenR$vectors
    Axis <- Re(ec[, abs(ev - 1) < 1e-10])
    if (abs(Axis[2])/max(svd(Axis)$d) > 1e-8) { naxis = xprod(Axis,c(1,0,0)) } ### Warning: not sure that norm(Axis) in Python corresponds to max(svd(Axis)$d) in R ###
    else { naxis = xprod(Axis,c(0,0,1)) }
    tn = R %*% naxis # i.e. matrix multiplication
    ca = tn %*% naxis
    sa = xprod(tn, naxis) %*% Axis
    a = atan2(sa, ca)
    saxis = Axis*a
    vm = matrix(c(0, saxis[3], -saxis[2], -saxis[3], 0, saxis[1], saxis[2], -saxis[1], 0), ncol = 3, byrow = T)
    return (list(vmoverDT = vm/DT, saxis = saxis))
    }
  else
    {
    stop("Cannot extract vorticity from a rotation in dimension greater than 3")
    }
  }


transformationToTensor <- function (t.M, DT, at_start = T, exp_correction = T)
  {
  if (nrow(t.M) != ncol(t.M)) { stop("The transformation is not an endomorphism") }
  if (nrow(t.M) > 3) { stop("This function doesn't work in more than 3 dimensions") }
  c(S, R) := polarDecomposition(t.M, at_start)
  R <- rotationToVorticity(R, DT)[[1]]
  S <- S - diag(nrow(S))
  S <- S / DT
  if (exp_correction) { S <- linearToExponential_strain(S, DT) }
  return (GT = S + R)
  }


growthTensor <- function (P, Q, DT, at_start = T, exp_correction = T, want_transform = F)
  {
  # Growth tensor transforming points P into points Q with DT.
  #
  #  The growth tensor will be the best one using least square criteria.

  t.M <- fitmat(P, Q)
  if (t.M[1] == "None") { return ("None") }
  GT <- transformationToTensor(t.M, DT, at_start, exp_correction)
  if (want_transform) { return (list(GT = GT, t.M = t.M)) }
  return (GT = GT)
  }


growthParams <- function (P, Q, DT, at_start = T, exp_correction = T)
  {
  # Return the growth parameters corresponding to the transformation of points
  # P into points Q with DT.
  #
  # :Parameters:
  #     P : ndarray(N*2)
  #         List of points at time t
  #     Q : ndarray(N*2)
  #         List of points at time t+dt
  #     DT : float
  #         Time between points P and Q
  #     exp_correction : bool
  #         If True, the result is corrected for exponential growth instead of linear
  #
  # :returns: The growth parameters as (kmaj, kmin, theta, phi) if 2D and
  #     (kmaj, kmed, kmin, theta_maj_xy, theta_maj_z, theta_med, psi_x, psi_y, psi_z) if 3D
  # :returntype: (float,)*4|(float,)*9

  GT <- growthTensor(P, Q, DT, at_start, exp_correction)
  if (GT[1] == "None")
    {
    return ("None")
    }
  else if (length(GT) == 1)
    {
    k <- GT[1]
    if (exp_correction) { return (linearToExponential_growth(k, DT)) }
    else { return (k) }
    }
  else if (nrow(GT) %in% c(2,3))
    {
    return (tensorToParams(GT))
    }
  else
    {
    stop("Cannot handle growth in more than 3D.")
    }
  }


tensorToParams <- function (tensor)
  {
  # :returns: The growth parameters as (kmaj, kmin, theta, phi) if 2D and
  #  (kmaj, kmed, kmin, theta_maj_xy, theta_maj_z, theta_med, psi_x, psi_y, psi_z) if 3D
  #  :returntype: (float,)*4|(float,)*9

  if (nrow(tensor) != ncol(tensor)) { stop("Tensor must be square.") }
  if (!nrow(tensor) %in% c(2, 3)) { stop("This function can only compute parameters for 2D and 3D tensors.") }

  # First, extract the symmetric and antisymmetric parts
  t.sym <- (tensor + t(tensor)) / 2
  t.anti <-(tensor - t(tensor)) / 2
  eigenTsym <- eigen(t.sym, symmetric = T) # returns eigenvalues in descending order
  values <- eigenTsym$values
  vectors <- eigenTsym$vectors
  # Eigenvalues and vectors should be sorted in the ascending order of the ABSOLUTE eigenvalues
  ascendingOrder <- order(abs(values))
  values <- values[ascendingOrder]
  vectors <- vectors[, ascendingOrder]

  if (nrow(tensor) == 2)
    {
    kmaj = values[2]
    kmin = values[1]
    theta = atan2(vectors[2,2], vectors[1,2])
    # Wrap over pi
    if (theta < pi/2) { theta <- theta + pi }
    else if (theta > pi/2) { theta <- theta - pi }
    phi <- t.anti[1,2]
    #if (round(abs(phi), 6) == round(pi, 6)) { phi <- 0 }
    return (list(kmaj = kmaj, kmin = kmin, theta = theta, phi = phi))
    }

  else   ### Warning: not tested ###
    {
    kmaj <- values[3]
    kmed <- values[2]
    kmin <- values[1]
    theta_maj_xy <- atan2(vectors[2,3], vectors[1,3])
    theta_maj_xyz <- atan2(vectors[3,3], max(svd(vectors[1:2,3])$d))
    theta_med <- atan2(max(svd(xprod(vectors[,2], vectors[,3]))$d), vectors[,2] %*% vectors[,3])
    # t.anti[i,j] = vorticity over vector i x j
    psi <- c(t.anti[2,3], -t.anti[1,3], t.anti[1,2])
    return (list(kmaj = kmaj, kmed = kmed, kmin = kmin, theta_maj_xy = theta_maj_xy, theta_maj_xyz = theta_maj_xyz, theta_med = theta_med, psi = psi))
    }
  }


paramsToTensor <- function(params)
  {
  # Convert growth parameters into growth tensor for 2D or 3D growth.

  if (length(params) == 4)
    {
    vx <- cos(params$theta)
    vy <- sin(params$theta)
    V.mat <- matrix(c(vx, vy, -vy, vx), nrow = 2, ncol = 2, byrow = T)
    D.mat <- diag(c(params$kmaj, params$kmin))
    Phi.mat <- matrix(c(0, params$phi, -params$phi, 0), nrow = 2, ncol = 2, byrow = T)
    return ((solve(V.mat) %*% D.mat %*% V.mat) + Phi.mat)
    }
  else if (length(params) == 9)
    {
    stop("The paramsToTensor() function is not yet implemented for 3D tensors.")
    }
  else
    {
    stop(paste("The paramsToTensor() function takes 4 or 9 arguments (", length(params), " given)", sep = ""))
    }
  }


################################################################################
#  R translation of Pierre Barbier de Reuille's growth_computation_methods.py  #
#        (partial, only for BackwardDenseMethod and ForwardDenseMethod)        #
################################################################################


length_polyline <- function (w)
  {
  x <- diff(w[, 1])
  y <- diff(w[, 2])
  return (sum(sqrt(x^2 + y^2)))
  }

length_segment <- function (S, Data)
  {
  total_length <- 0
  lengths <- 0
  for (p.iter in 1:(length(S)-1))
    {
    p1 <- S[p.iter]
    p2 <- S[p.iter+1]
    ################### w <- data.walls[p1,p2]
    w <- matrix(nrow = 0, ncol = 2) ###################
    w <- rbind(Data[row.names(Data) == p1, ], w)
    w <- rbind(w, Data[row.names(Data) == p2, ])
    l <- length_polyline(w)
    total_length <- total_length + l
    lengths <- c(lengths, total_length)
    }
  return (lengths)
  }

subfn_align <- function (Length, S, Ratios, Data, all_pos, len_s2) # (Originally a subfunction of 'align_segments()')
  {
  align <- matrix(ncol = 3, nrow = (length(all_pos)-1))
  Ratios <- unique(Ratios)
  p1 <- 1 # Position in s1
  p2 <- 2 # Position in s1
  pos <- matrix(Data[row.names(Data) == S[p1], ], ncol = 2)
  cur <- S[p1]
  Next <- NA
  ################### w <- data.walls[S[p1],S[p2]]
  w <- matrix(nrow = 0, ncol = 2) ###################

  for (j in 1:(length(all_pos)-1))
    {
    r1 <- all_pos[j]
    r2 <- all_pos[j+1]
    w1 <- pos

    if (r2 %in% Ratios) # If the next point is in the current wall
      {
      w1 <- rbind(w1, w)
      pos <- matrix(Data[row.names(Data) == S[p2], ], ncol = 2)
      p1 <- p1+1
      p2 <- p2+1
      if (p2 < (length(S)+1))
        {
        ################### w <- data.walls[S[p1],S[p2]]
        w <- matrix(nrow = 0, ncol = 2) ###################
        Next <- S[p1]
        }
      }

    else # Otherwise, find where it stops
      {
      l <- (r2-r1)*Length
      acc <- 0
      while (nrow(w) > 0)
        {
        p <- w[1, ]
        vec <- p - w1[nrow(w1), ]
        dl <- sqrt(vec[,1]^2 + vec[,2]^2)
        if (acc + dl > l) # If we go past the next stop
          {
          pos <- w1[nrow(w1), ] + vec*((l-acc)/dl)
          break
          }
        acc <- acc + dl
        w1 <- rbind(w1, p)
        w <- w[-1, ]
        }
      if (nrow(w) == 0) # If the end of the wall has been reached
        {
        p <- matrix(Data[row.names(Data) == S[p2], ], ncol = 2)
        vec <- p - w1[nrow(w1), ]
        dl <- sqrt(vec[,1]^2 + vec[,2]^2)
        pos <- w1[nrow(w1), ] + vec*((l-acc)/dl)
        }
      }

    align[j, ] <- c(w1, cur)
    l <- (r2-r1)*len_s2
    cur <- Next
    Next <- NA
    }

  return (align)
  }

align_segments <- function (s1, s2, data1, data2)
  {
  # Compute the alignment of segments s1 and s2,
  # such that the first and last elements of s1 and s2 are the same, but nothing else.
  #
  # :return_type: list of (list of QPointF, int)
  # :returns: List of wall parts such that the first point is the vertex.
  #           The integer is the id of the point (if it corresponds to one).

  # First, compute ratios
  lengths_s1 <- length_segment(s1, data1)
  lengths_s2 <- length_segment(s2, data2)
  ratios_s2 <- ratios_s1 <- vector()
  for (l in lengths_s1) { ratios_s1 <- c(ratios_s1, l/lengths_s1[length(lengths_s1)]) }
  for (l in lengths_s2) { ratios_s2 <- c(ratios_s2, l/lengths_s2[length(lengths_s2)]) }
  len_s1 <- lengths_s1[length(lengths_s1)]
  len_s2 <- lengths_s2[length(lengths_s2)]
  all_pos <- sort(unique(c(ratios_s1, ratios_s2)))
  align_s1 <- subfn_align(len_s1, s1, ratios_s1, data1, all_pos, len_s2)
  align_s2 <- subfn_align(len_s2, s2, ratios_s2, data2, all_pos, len_s2)
  return (list(align_s1 = align_s1, align_s2 = align_s2))
  }

discretize_segment <- function (seg, n, l)
  {
  if (n == 0) # two points have been placed exactly at the same position (Image1)
    {
    result <- matrix(seg[1,], nrow = 1, ncol = 2, byrow = T)
    rownames(result) <- "pos"
    }  
  else if (l == 0) # two points have been placed exactly at the same position (Image2)
    {
    result <- matrix(seg[1,], nrow = n, ncol = 2, byrow = T)
    rownames(result) <- rep("pos", n)
    }
  else
    {
    dl <- l/n
    result <- matrix(nrow = 0, ncol = 2)
    idx <- 1
    pos <- seg[idx, ]
    vec <- seg[idx+1, ] - pos
    vec_size <- sqrt(vec[1]^2 + vec[2]^2)
    vec <- vec/vec_size
    shift <- 0
    for (j in 1:n)
      {
      result <- rbind(result, pos)
      if (j == n) { break }
      needed_dl <- dl
      while (shift + needed_dl > vec_size)
        {
        idx <- idx+1
        needed_dl <- needed_dl - (vec_size - shift)
        pos <- seg[idx, ]
        vec <- seg[idx+1, ] - seg[idx, ]
        vec_size <- sqrt(vec[1]^2 + vec[2]^2)
        vec <- vec/vec_size
        shift <- 0
        }
      pos <- pos + vec*needed_dl
      shift <- shift + needed_dl
      }
    }
  return (result)
  }

alignCells <- function (cell, pts, new_pts, img_data, next_img_data, nb_points, Image1, Image2)
  {
  # First, find a common vertex between the cells
  error = T
  for (common in pts)
    {
    if (common %in% new_pts)
      {
      idx <- which(pts == common)
      idx1 <-  which(new_pts == common)
      pts <- c(pts[idx:length(pts)], pts[0:(idx-1)])
      new_pts <- c(new_pts[idx1:length(new_pts)], new_pts[0:(idx1-1)])
      error = F
      break
      }
    }
  if (error)
    {
    stop(paste("Error,", cell, "have no common point between times", Image1, "and", Image2))
    }

  # Then, align the cells. i.e. add missing points
  aligned_pts <- aligned_new_pts <- matrix(nrow = 0, ncol = 3)
  i1 <- j1 <- 1
  for (i2 in 2:length(pts))
    {
    pid <- pts[i2]
    if (pid %in% new_pts)
      {
      j2 <- which(new_pts == pid)
      if (j2 < j1)
        {
        cat(paste("Error:", cell, "is inconsistent between times", Image1, "and", Image2))
        aligned_pts <- aligned_new_pts <- vector()
        i1 <- j1 <- 1
        break
        }
      c(seg, new_seg) := align_segments(pts[i1:i2], new_pts[j1:j2], img_data, next_img_data)
      aligned_pts <- rbind(aligned_pts, seg)
      aligned_new_pts <- rbind(aligned_new_pts, new_seg)
      i1 <- i2
      j1 <- j2
      }
    }

  # Add the missing segment
  if (i1 != 1)
    {
    c(seg, new_seg) := align_segments(c(pts[i1:length(pts)], pts[1]), c(new_pts[j1:length(new_pts)], new_pts[1]), img_data, next_img_data)
    aligned_pts <- rbind(aligned_pts, seg)
    aligned_new_pts <- rbind(aligned_new_pts, new_seg)
    }
  if (nrow(aligned_pts) == 0) { return () }

  # Next, for the cell, start by resampling them
  # Compute total perimeter of first cell
  l1 <- aligned_pts[, 1:2]
  l1 <- rbind(l1, l1[1,])
  len_c1 <- length_polyline(l1)
  # Approximate the dl
  dl <- len_c1 / nb_points
  Ps <- matrix(nrow = 0, ncol = 2) # Resampled first cell
  Qs <- matrix(nrow = 0, ncol = 2) # Resampled second cell
  nb_seg <- nrow(aligned_pts)
  for (i in 1:nrow(aligned_pts))
    {
    seg1 <- aligned_pts[i, 1:2]
    seg2 <- aligned_new_pts[i, 1:2]
    seg1 <- rbind(seg1, aligned_pts[(i %% nb_seg) + 1, 1:2])
    seg2 <- rbind(seg2, aligned_new_pts[(i %% nb_seg) + 1, 1:2])
    l1 <- length_polyline(seg1)
    l2 <- length_polyline(seg2)
    # Number of points on the current segment
    n <- ceiling(l1/dl)
    # Real dl for first and second cells
    Ps <- rbind(Ps, discretize_segment(seg1, n, l1))
    Qs <- rbind(Qs, discretize_segment(seg2, n, l2))
    }

  return (list(aligned_pts = aligned_pts, aligned_new_pts = aligned_new_pts, Ps = Ps, Qs = Qs))
  }



################################################################################
#          R translation of Richard Kennaway's color scale for GFtbox          #
################################################################################

  
rainbowMap <- function (Range, zerowhite = F, n.colors)
  {
  negpart <- matrix(c(0.9, 0, 1,     # purple
                      0.45, 0, 0.5), # dark purple
                      nrow = 2, byrow = T)

  pospart <- matrix(c(0, 0, 1,     # blue
                      0, 0.7, 1,   # bluish cyan
                      0, 1, 1,     # cyan
                      0, 1, 0.7,   # greenish cyan
                      0, 1, 0,     # green
                      0.75, 1, 0,  # chartreuse
                      1, 1, 0,     # yellow
                      1, 0.875, 0, # yellow-orange
                      1, 0.5, 0,   # orange
                      1, 0, 0,     # red
                      2/3, 0, 0),  # darker red
                      nrow = 11, byrow = T)
  
  if (zerowhite) { pospart <- rbind(c(1, 1, 1), pospart) } # fades to white at zero

  negpart <- rbind(pospart[1,], negpart)
    
  if (missing(Range)) { Range <- c(-1, 1) }
  #Range <- extendToZero(Range) ## Not used since it must be done only if missing(zlim); set.zlim() does the job
  Colors <- posnegMap(Range, negpart, pospart, n.colors)
  return (list(Colors = rgb(Colors), Range = Range))
  }
  
extendToZero <- function (Range)
  {
  if (length(Range) != 0)
    {
    if (length(Range) == 1)
      {
      ifelse (Range <= 0, Range <- c(Range, 0), Range <- c(0, Range))
      }
    else
      {
      if (Range[1] >= 0) { Range[1] <- 0 }
      if (Range[2] <= 0) { Range[2] <- 0 }
      }
    }
  return (Range)
  }


posnegMap <- function (Range, negpart, pospart, n.colors)
  {
  if (Range[1] == Range[2])
    {
    Colors <- matrix(rep(1,6), nrow = 2, byrow = T)
    return (Colors)
    }

  nsteps <- n.colors
  negfrac <- Range[1]/(Range[1] - Range[2])
  posfrac <- Range[2]/(Range[2] - Range[1])
  nneg <- max(c(0, ceiling(negfrac*nsteps)))
  npos <- max(c(0, ceiling(posfrac*nsteps)))
  maxfrac <- max(c(negfrac, posfrac))
  negfrac <- negfrac/maxfrac
  posfrac <- posfrac/maxfrac
  if (nneg == 0)
    {
    posmap <- makeCmap(pospart, npos+1, posfrac)
    Colors <- posmap
    }
  else if (npos == 0)
    {
    negmap <- makeCmap(negpart, nneg, negfrac)
    Colors <- negmap[nrow(negmap):1, ]
    }
  else
    {
    negmap <- makeCmap(negpart, nneg, negfrac)
    posmap <- makeCmap(pospart, npos, posfrac)
    Colors <- rbind(negmap[nrow(negmap):1, ], posmap[2:nrow(posmap), ])
    }

  return (Colors)
  }

makeCmap <- function (Col, nsteps, frac)
  {
  if (nsteps <= 0)
    {
    cmap <- matrix(0, 0, 3)
    return (cmap)
    }

  cmap <- matrix(0, nsteps, 3)
  ncsteps <- nrow(Col)-1
  a <- (0:(nsteps))*frac
  b <- a*(ncsteps/nsteps)
  c <- floor(b)
  d <- b-c
  e <- 1-d
  for (i in 1:nsteps)
    {
    ci <- c[i]
    cmap[i, ] <- Col[ci+1,]*e[i] + Col[ci+2,]*d[i]
    }

  return (cmap)
  }
  

colorMap.Residual <- function (Range, zerowhite = T, n.colors)
  {
  # Blue-white-red colour scale for residual data

  negpart <- matrix(c(0, 1, 1,    # cyan
                      0, 0.7, 1,   # bluish cyan
                      0, 0, 1),     # blue
                      nrow = 3, byrow = T)

  pospart <- matrix(c(1, 0.5, 0,   # orange
                      1, 0, 0,     # red
                      2/3, 0, 0),  # darker red
                      nrow = 3, byrow = T)

  if (zerowhite) { pospart <- rbind(c(1, 1, 1), pospart) } # fades to white at zero

  negpart <- rbind(pospart[1,], negpart)

  if (missing(Range)) { Range <- c(-1, 1) }
  #Range <- extendToZero(Range) ## Not used since it must be done only if missing(zlim); set.zlim() does the job
  Colors <- posnegMap(Range, negpart, pospart, n.colors)
  return (list(Colors = rgb(Colors), Range = Range))
  }

  
  
################################################################################
#                        Leaf shape and cell coordinates                       #
################################################################################


fn.rotation <- function (matXY, Angle)
  {
  # This function returns the rotated coordinates of an input shape.
  # Xs should be in the first column, and Ys should be in the second column.
  # 'Angle' should be given in radians.
  
  xRot <- matXY[,1] * cos(Angle) + matXY[,2] * sin(Angle)
  yRot <- -matXY[,1] * sin(Angle) + matXY[,2] * cos(Angle)
  return (cbind(xRot, yRot))
  }


get.meanAngle <- function (InfoVertices)
  {
  # This function returns the mean angle (in degrees) by which all images will be rotated
  # so that each leaf will be oriented upward. Requires digitising leaf outline using ImageJ.
  # It takes into account that the leaf will be also rotated by AngleShift
  # as computed in PointTracker.
  
  angle_leaf <- vector()
  for (img in InfoVertices$Images)
    {
    LeafShape <- read.table(paste("Leaf Shape/", substr(img, 1, nchar(img) - 3), "txt", sep = ""))
    centro <- get.Centroid(LeafShape)
    if (centro[1] == mean(LeafShape[1:2, 1]))
      {
      ifelse (min(LeafShape[, 1]) %in% LeafShape[1:2, 1], leaf_orientation <- 90, leaf_orientation <- -90)
      }
    else
      {
      leaf_orientation <- -atan2(centro[2] - mean(LeafShape[1:2, 2]), centro[1] - mean(LeafShape[1:2, 1])) * 180/pi # if close to 90 leaf is vertical, if 0 leaf is horizontal, ...
      }
    AngleShift <- -InfoVertices$AngleShift[InfoVertices$Images == img]
    angle_leaf <- c(angle_leaf, 90 - (leaf_orientation + AngleShift))
    }
  return (mean(angle_leaf))
  }
  
  
get.leafShape <- function (Image, meanAngle, InfoVertices)
  {
  # This function retrieves leaf outline from ImageJ
  # and standardize coordinates with those of PointTracker.
  # As in the ImageJ macro 'Draw_leaf_outlines.txt',
  # the image should not have been scaled in microns
  # and should not have been rotated.
  # The two first points should form the petiole base.

  LeafShape <- read.table(paste("Leaf Shape/", substr(Image, 1, nchar(Image) - 3), "txt", sep = ""))
  scaleX <- InfoVertices$ScalingX[InfoVertices$Images == Image]
  scaleY <- InfoVertices$ScalingY[InfoVertices$Images == Image]
  AngleShift <- -InfoVertices$AngleShift[InfoVertices$Images == Image]
  LeafShape$V1 <- LeafShape$V1 * scaleX * 1e6
  LeafShape$V2 <- LeafShape$V2 * scaleY * 1e6
  Angle_Rotation <- (AngleShift + meanAngle) * pi/180
  LeafShape[, c("V1", "V2")] <- fn.rotation(LeafShape[, c("V1", "V2")], Angle_Rotation)
  LeafShape$X <- LeafShape$V1 - mean(LeafShape$V1[1:2])
  #LeafShape$Y <- mean(LeafShape$V2[1:2]) - LeafShape$V2
  LeafShape$Y <- min(LeafShape$V2[1:2]) - LeafShape$V2
  LeafShape$Y[1:2] <- 0 # allows straight petiole
  return (LeafShape)
  }


get.leafSpline <- function (LeafShape)
  {
  # This function returns a spline from the leaf outline,
  # with a straight (horizontal) petiole base
  
  current.device <- dev.cur()
  windows()
  plot(0)
  LeafSpline <- xspline(LeafShape$X, LeafShape$Y, s = c(0, 0, rep(1, nrow(LeafShape) - 2)), open = F, draw = F)
  dev.off()
  if (current.device != 1) { dev.set(current.device) }
  return (LeafSpline)
  }


process.Vertices <- function (Image, LeafShape, meanAngle, InfoVertices, Vertices, InputVertices)
  {
  # This function retrieves the vertices from PointTracker at a given image
  # and transform the coordinates to match leaf outline.
  # Vertices can be given as 'InputVertices' (eg as an individual cell 'Data',
  # i.e. a densified shape from 'cell.coord.for.alignCells()',
  # but 'alignToPLB()' won't work on the output vertices
  # because the latter function requires all the cells to compute the offset).
  # 'meanAngle' shall be given in degrees.

  scaleX <- InfoVertices$ScalingX[InfoVertices$Images == Image]
  scaleY <- InfoVertices$ScalingY[InfoVertices$Images == Image]
  shiftX <- InfoVertices$ShiftX[InfoVertices$Images == Image]
  shiftY <- InfoVertices$ShiftY[InfoVertices$Images == Image]

  if (missing(InputVertices))
    {
    VerticesX <- as.numeric(Vertices[, which(names(Vertices) == Image)[1]])
    VerticesY <- as.numeric(Vertices[, which(names(Vertices) == Image)[2]])
    }
  else
    {
    VerticesX <- InputVertices[,1]
    VerticesY <- InputVertices[,2]
    }
  VerticesX <- (VerticesX - shiftX * scaleX) * 1e6
  VerticesY <- (VerticesY - shiftY * scaleY) * 1e6

  rotated.vert <- fn.rotation(cbind(VerticesX, VerticesY), meanAngle * pi/180)
  VerticesX <- rotated.vert[,1]
  VerticesY <- rotated.vert[,2]

  VerticesX <- VerticesX - mean(LeafShape$V1[1:2])
  #VerticesY <- mean(LeafShape$V2[1:2]) - VerticesY
  VerticesY <- min(LeafShape$V2[1:2]) - VerticesY

  return (list(VerticesX = VerticesX, VerticesY = VerticesY))
  }


process.Shapes <- function (Image, before, LeafShape, meanAngle, InfoVertices, Shapes)
  {
  # This function retrieves the cell shapes from PointTracker at a given image
  # and transform the coordinates to match leaf shape.
  # Set before to TRUE if 'Image' is the initial step, and to FALSE if 'Image' is the final step.
  # 'meanAngle' shall be given in degrees.
   
  scaleX <- InfoVertices$ScalingX[InfoVertices$Images == Image]
  scaleY <- InfoVertices$ScalingY[InfoVertices$Images == Image]
  shiftX <- InfoVertices$ShiftX[InfoVertices$Images == Image]
  shiftY <- InfoVertices$ShiftY[InfoVertices$Images == Image]

  ifelse (before, row.offset <- 1, row.offset <- 0)

  ShapesX <- Shapes[1:(nrow(Shapes)/2)*2 - row.offset, c(1, seq(3, ncol(Shapes)-1, 2))]
  ShapesY <- Shapes[1:(nrow(Shapes)/2)*2 - row.offset, c(1, seq(3+1, ncol(Shapes), 2))]

  ShapesX[, 2:ncol(ShapesX)] <- (ShapesX[, 2:ncol(ShapesX)] - shiftX * scaleX) * 1e6
  ShapesY[, 2:ncol(ShapesY)] <- (ShapesY[, 2:ncol(ShapesY)] - shiftY * scaleY) * 1e6

  xRot <- ShapesX[, 2:ncol(ShapesX)] * cos(meanAngle * pi/180) + ShapesY[, 2:ncol(ShapesY)] * sin(meanAngle * pi/180)
  yRot <- -ShapesX[, 2:ncol(ShapesX)] * sin(meanAngle * pi/180) + ShapesY[, 2:ncol(ShapesY)] * cos(meanAngle * pi/180)
  ShapesX[, 2:ncol(ShapesX)] <- xRot
  ShapesY[, 2:ncol(ShapesY)] <- yRot

  ShapesX[, 2:ncol(ShapesX)] <- ShapesX[, 2:ncol(ShapesX)] - mean(LeafShape$V1[1:2])
  #ShapesY[, 2:ncol(ShapesY)] <- mean(LeafShape$V2[1:2]) - ShapesY[, 2:ncol(ShapesY)]
  ShapesY[, 2:ncol(ShapesY)] <- min(LeafShape$V2[1:2]) - ShapesY[, 2:ncol(ShapesY)]
  
  return (list(ShapesX = ShapesX, ShapesY = ShapesY))
  }


alignToPLB <- function (VerticesY, cellAtPetioleLaminaBoundary, Cells, Vertices)
  {
  # This function returns the Y-offset to align the vertices and the leaf shape
  # to the petiole-lamina boundary.
  # 'cell' should be provided as an integer being the cell at the petiole-lamina boundary.
  # If no 'cell' is provided, the lowest tracked vertex is used.
  
  if (missing(cellAtPetioleLaminaBoundary))
    {
    Offset <- min(VerticesY, na.rm = T)
    }
  else
    {
    if (length(cellAtPetioleLaminaBoundary) != 1 | as.integer(cellAtPetioleLaminaBoundary) != cellAtPetioleLaminaBoundary | !paste("Cell", cellAtPetioleLaminaBoundary) %in% Cells$Cell)
      {
      stop("The parameter 'cellAtPetioleLaminaBoundary' should be an integer corresponding to an existing cell.")
      }
    ifelse (is.null(Divisions), cellsAtPetioleLaminaBoundary <- cellAtPetioleLaminaBoundary,
                                cellsAtPetioleLaminaBoundary <- find.Lineage(cellAtPetioleLaminaBoundary, Divisions))
    VerticesOfCellsAtPetioleLaminaBoundary <- paste("Point", unique(na.omit(as.numeric(as.matrix(Cells[Cells$Cell %in% cellsAtPetioleLaminaBoundary, -1])))))
    Offset <- min(VerticesY[Vertices$Point %in% VerticesOfCellsAtPetioleLaminaBoundary], na.rm = T)
    }
  return (Offset)
  }


cell.coord <- function (vert, VerticesX, VerticesY, Vertices)
  {
  # This function returns the coordinates of the current cell
  # given a set of vertices 'vert'.

  matXY <- matrix(c(VerticesX[Vertices$Point %in% paste("Point", vert)], VerticesY[Vertices$Point %in% paste("Point", vert)]), ncol = 2)
  matXY <- matXY[rank(vert), ]
  if (length(which(is.na(matXY[,1]))) > 0) { matXY <- matXY[-which(is.na(matXY[,1])), ] }
  return (matXY)
  }


cell.coord.for.alignCells <- function (Image, vert, Vertices)
  {
  # This function returns the vertices and coordinates required by 'alignCells()'
  # given a set of vertices 'vert'.

  # Coordinates of the vertices at Image (ascending order, including NAs)
  Data <- as.matrix(Vertices[Vertices$Point %in% paste("Point", vert), which(names(Vertices) == Image)])
  # ID of the vertices at Image (i.e. 'vert' excluding NAs of 'Data')
  ifelse (length(which(is.na(Data[,1]))) > 0, pts <- vert[-order(vert)[which(is.na(Data[,1]))]], pts <- vert)
  # Rearranges 'Data' for 'vert'
  Data <- Data[rank(vert), ]
  # Rearranges 'Data' for 'pts'
  if (length(which(is.na(Data[,1]))) > 0) { Data <- Data[-which(is.na(Data[,1])), ] }
  row.names(Data) <- pts
  return (list(Data = Data, pts = pts))
  }


get.Centroid <- function (matXY)
  {
  # This function returns the (X,Y) coordinates of the centroid of a polygon
  # using the 'gCentroid' function of the package 'rgeos'.

  poly.string <- "POLYGON(("
  for (v in 1:nrow(matXY)) { poly.string <- paste(poly.string, matXY[v,1], " ", matXY[v,2], ",", sep = "") }
  poly.string <- paste(poly.string, matXY[1,1], " ", matXY[1,2], "))", sep = "") # close the polygon
  poly.wkt <- readWKT(poly.string)
  centroid <- gCentroid(poly.wkt)
  return (centroid@coords)
  }


scale.Objects <- function (Image1, Image2, before, meanAngle, leafShape, alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary, Shapes, Cells, Vertices, InfoVertices)
  {
  # This function prepares the objects to be plotted in a coordinate system
  # where axes are graduated in microns and the origin is the intersection of the midvein with:
  #   - if 'alignToPetioleLaminaBoundary' is TRUE, the y-coordinate of lowest vertex
  #     of the cell corresponding to the petiole-lamina boundary;
  #   - if 'alignToPetioleLaminaBoundary' is FALSE, the y-coordinate of the base
  #     of the leaf outline (i.e. the petiole).
  # The y-axis is parallel to the midvein if 'meanAngle' is provided accordingly
  # (i.e. if computed using 'get.meanAngle()' which works if the leaf outline is provided).

  ifelse (before, Image <- Image1, Image <- Image2)
  
  if (leafShape)
    {
    LeafShape <- get.leafShape(Image, meanAngle, InfoVertices)
    LeafSpline <- get.leafSpline(LeafShape)
    }
  else
    {
    LeafShape <- LeafShape <- data.frame(V1 = c(0, 0), V2 = c(0, 0), X = c(0, 0), Y = c(0, 0))
    LeafSpline <- LeafSpline <- list(x = 0, y = 0)
    }
    
  c(VerticesX, VerticesY) := process.Vertices(Image, LeafShape, meanAngle, InfoVertices, Vertices)
  c(ShapesX, ShapesY) := process.Shapes(Image, before, LeafShape, meanAngle, InfoVertices, Shapes)

  if (alignToPetioleLaminaBoundary)
    {
    Offset <- alignToPLB(VerticesY, cellAtPetioleLaminaBoundary, Cells, Vertices)
    VerticesY <- VerticesY - Offset
    LeafShape$Y <- LeafShape$Y - Offset
    LeafSpline$y <- LeafSpline$y - Offset
    ShapesY[, 2:ncol(ShapesY)] <- ShapesY[, 2:ncol(ShapesY)] - Offset
    }
  
  return (list(VerticesX = VerticesX, VerticesY = VerticesY, ShapesX = ShapesX, ShapesY = ShapesY, LeafShape = LeafShape, LeafSpline = LeafSpline))
  }
   


################################################################################
#                                   Lineages                                   #
################################################################################


find.Cells <- function (Image, Cells, Divisions, Vertices, InfoVertices)
  {
  # This function finds the current cells at the image considered.
  
  # --> edit 2016/09/22 and 2017/02/24 #
  InfoVertices$Time <- round(InfoVertices$Time, 4)
  if (!is.null(Divisions))
    {
    Divisions$TimeHours1 <- round(Divisions$TimeHours1, 4)
    Divisions$TimeHours2 <- round(Divisions$TimeHours2, 4)
    }
  # edit 2016/09/22 and 2017/02/24 <-- #

  Vertices.Image <- as.character(Vertices$Point[!is.na(Vertices[, which(names(Vertices) == Image)[1]])])
  Vertices.Image <- as.numeric(substr(Vertices.Image, 7, nchar(Vertices.Image)))
  levels.Cells <- vector()
  for (cell in 1:nrow(Cells))
    {
    if (length(which(Cells[cell, -1] %in% Vertices.Image)) > 2) # At least 3 points to define a cell
      {
      if (is.null(Divisions))
        {
        levels.Cells <- c(levels.Cells, as.character(Cells$Cell[cell])) 
        }
      else
        {        
        check.Daughter <- check.Mother <- T
        if (as.character(Cells$Cell[cell]) %in% Divisions$Cell) # If the cell has divided...
          {
          if (Divisions$TimeHours2[Divisions$Cell == as.character(Cells$Cell[cell])] <= InfoVertices$Time[InfoVertices$Image == Image]) #... was it before the image considered?
            {
            check.Mother <- F
            }
          }
        if (as.character(Cells$Cell[cell]) %in% paste("Cell", c(Divisions[, "Daughter cell 1"], Divisions[, "Daughter cell 2"]))) # If the cell is a daughter...
          {
          if (rep(Divisions$TimeHours2, 2)[paste("Cell", c(Divisions[, "Daughter cell 1"], Divisions[, "Daughter cell 2"])) == as.character(Cells$Cell[cell])] > InfoVertices$Time[InfoVertices$Image == Image]) #... was it born after the image considered?
            {
            check.Daughter <- F
            }
          }
        if (check.Daughter & check.Mother)
          {
          levels.Cells <- c(levels.Cells, as.character(Cells$Cell[cell]))
          }
        }
      }
    }
  return (levels.Cells)
  }


find.OldestAncestor <- function (cell, time1, Divisions)
  {
  # This function returns the direct mother of a cell and its oldest ancestor
  # at 'time1' (which should be in the same unit as in PointTracker, i.e. hours). 

  # --> edit 2017/02/26 #
  if (is.null(Divisions))
    {
    oldest.ancestor <- cell
    direct.mother <- vector()
    n.div <- 0
    }
  # edit 2017/02/26 <-- #

  else
    {
    # --> edit 2016/09/22 #
    time1 <- round(time1, 4)
    Divisions$TimeHours1 <- round(Divisions$TimeHours1, 4)
    Divisions$TimeHours2 <- round(Divisions$TimeHours2, 4)
    # edit 2016/09/22 <-- #

    grand.mother <- direct.mother <- as.character(rep(Divisions$Cell, 2)[which(paste("Cell", c(Divisions[, "Daughter cell 1"], Divisions[, "Daughter cell 2"])) == cell)])
    oldest.ancestor <- cell
    n.div <- 0 # edit 2015/12/15

    while (length(grand.mother) == 1)
      {
      grand.mother <- as.character(rep(Divisions$Cell, 2)[which(paste("Cell", c(Divisions[, "Daughter cell 1"], Divisions[, "Daughter cell 2"])) == oldest.ancestor)])
      if (length(grand.mother) == 1)
        {
        if (Divisions$TimeHours2[Divisions$Cell == grand.mother] > time1) # was the mother born after time1?
          {
          oldest.ancestor <- grand.mother
          n.div <- n.div + 1 # edit 2015/12/15
          }
        else
          {
          grand.mother <- character(0)
          }
        }
      }
    }
  return (list(oldest.ancestor = oldest.ancestor, direct.mother = direct.mother, n.div = n.div)) # edit 2015/12/15
  }
 

find.Lineage <- function (target.cell, Divisions)
  {
  # This function finds the lineage of a target cell.
  # The oldest ancestor is found, and all descendants of this ancestor are returned.
  # The target cell should be given as an integer.
  
  c(mother, direct.mother) := find.OldestAncestor(paste("Cell", target.cell), 0, Divisions)
  
  lineage <- mother
  new.mother <- vector()
  check.lineage <- F
  while (check.lineage == F)
    {
    for (m in mother)
      {
      D1 <- Divisions[Divisions$Cell == m, "Daughter cell 1"]
      D2 <- Divisions[Divisions$Cell == m, "Daughter cell 2"]
      if (length(D1) == 1)  
        {
        lineage <- c(lineage, paste("Cell", D1), paste("Cell", D2))
        new.mother <- c(new.mother, paste("Cell", D1), paste("Cell", D2))
        }
      else if (m == mother[length(mother)] & length(new.mother) == 0)
        {
        check.lineage <- T
        }
      }
    mother <- new.mother
    new.mother <- vector()
    }
  
  return (lineage)
  }


# --> edit 2015/12/15 #
find.MaxNumberOfDivisions <- function (cell, time1, time2, Divisions, before)
  {
  # This function returns the maximal number of divisions for a cell between 'time1' and 'time2'.
  # If 'before' is FALSE, then the maximal number of divisions equals the number of divisions given by 'find.OldestAncestor()'.
  # If 'before' is TRUE, then the maximal number of divisions equals is unique to the oldest ancestor at 'time1'.

  # --> edit 2016/09/22 #
  time1 <- round(time1, 4)
  time2 <- round(time2, 4)
  Divisions$TimeHours1 <- round(Divisions$TimeHours1, 4)
  Divisions$TimeHours2 <- round(Divisions$TimeHours2, 4)
  # edit 2016/09/22 <-- #

  if (!before)
    {
    c(mother, direct.mother, n.div) := find.OldestAncestor(cell, time1, Divisions[Divisions$TimeHours2 > time1 & Divisions$TimeHours2 <= time2, ])
    }

  else
    {
    n.div <- vector()
    lineage <- find.Lineage(as.numeric(substr(cell, 6, nchar(cell))), Divisions[Divisions$TimeHours2 > time1 & Divisions$TimeHours2 <= time2, ])
    for (l in lineage)
      {
      c(mother, direct.mother, n.div.temp) := find.OldestAncestor(l, time1, Divisions[Divisions$TimeHours2 > time1 & Divisions$TimeHours2 <= time2, ])
      n.div <- c(n.div, n.div.temp)
      }
    n.div <- max(n.div)
    }

  return (n.div)
  }
# edit 2015/12/15 <-- #
   
################################################################################
#                         Import data from PointTracker                        #
################################################################################


import.rawData <- function(File)
  {
  # This function imports the raw data from PointTracker.
  # The csv file should be in a subfolder of the working directory called 'Data'.
  #
  # Below is an example for the 'File' argument:
  #     File = "tracking.csv"

  # Import csv file
  dataImport <- read.table(paste("Data", File, sep = "/"), sep = ",", header = F, fill = T, blank.lines.skip = F, col.names = 1:1000, na.strings = "") # edit 2018/07/07
  activeCol <- vector()
  for (Col in 1:ncol(dataImport)) { activeCol <- c(activeCol, length(which(is.na(dataImport[Col])))) }
  if (length(which(activeCol == nrow(dataImport))) != 0) { dataImport <- dataImport[-which(activeCol == nrow(dataImport))] }

  # Image properties stored in the csv file
  Info <- dataImport[which(dataImport[1] == "Images") + 0:3, 1:length(which(!is.na(dataImport[(which(dataImport[1] == "Images") + 3), ])))]
  InfoVertices <- data.frame(Images = NULL, ShiftX = NULL, ShiftY = NULL, AngleShift = NULL, Time = NULL, ScalingX = NULL, ScalingY = NULL)
  for (imag in seq(2, ncol(Info), 2)) { InfoVertices <- rbind(InfoVertices, data.frame(Images = as.character(Info[1, imag]), ShiftX = as.numeric(as.character(Info[2, imag])), ShiftY = as.numeric(as.character(Info[2, imag+1])), AngleShift = as.numeric(as.character(Info[3, imag])), Time = as.numeric(as.character(Info[3, imag+1])), ScalingX = as.numeric(as.character(Info[4, imag])), ScalingY = as.numeric(as.character(Info[4, imag+1])))) }
  InfoVertices$Images <- as.factor(InfoVertices$Images)

  # Vertices coordinates
  Vertices <- dataImport[(which(dataImport[1] == "Images") + 4):(which(dataImport[1] == "Cells") - 1), 1:length(which(!is.na(dataImport[(which(dataImport[1] == "Images") + 4), ])))]
  names(Vertices)[1] <- "Point"
  for (Col in seq(2, ncol(Vertices), 2)) { names(Vertices)[c(Col, Col+1)] <- as.character(dataImport[(which(dataImport[1] == "Images") + 0), Col]) }
  Vertices[, 1] <- factor(Vertices[, 1])
  for (Col in 2:ncol(Vertices)) { Vertices[, Col] <- as.numeric(as.character(Vertices[, Col])) }
     
  # Cells to vertices
  Cells <- dataImport[(which(dataImport[1] == "Cells") + 1):(which(dataImport[1] == "Divisions") - 1), ]
  Cells[, 1] <- factor(Cells[, 1])
  for (Col in 2:ncol(Cells)) { Cells[, Col] <- as.numeric(as.character(Cells[, Col])) }
  activeCol <- vector()
  for (Col in 1:ncol(Cells)) { activeCol <- c(activeCol, length(which(is.na(Cells[Col])))) }
  if (length(which(activeCol == nrow(Cells))) != 0) { Cells <- Cells[-which(activeCol == nrow(Cells))] }
  names(Cells)[1] <- "Cell"
  names(Cells)[2:ncol(Cells)] <- paste("Vertex", 2:ncol(Cells) - 1, sep = "")
     
  # Divisions
  if (which(dataImport[1] == "Divisions") < nrow(dataImport))
    {
    Divisions <- dataImport[(which(dataImport[1] == "Divisions") + 1):nrow(dataImport), 1:6]
    names(Divisions)[1] <- "Cell"
    for (Col in 2:ncol(Divisions)) { Divisions[, Col] <- as.numeric(as.character(Divisions[, Col])) }
    for (Col in 2:ncol(Divisions)) { names(Divisions)[Col] <- as.character(dataImport[(which(dataImport[1] == "Divisions")), Col]) }
    for (t in as.numeric(as.character(levels(as.factor(Divisions$Time)))))
      {
      Divisions$TimeHours1[Divisions$Time == t] <- InfoVertices$Time[t]    # Time before cell division (hours)
      Divisions$TimeHours2[Divisions$Time == t] <- InfoVertices$Time[t+1]  # Time after cell division (hours)
      }
    Divisions$Cell <- factor(Divisions$Cell)
    }
  else
    {
    Divisions <- NULL
    }
  
  return (list(InfoVertices = InfoVertices, Vertices = Vertices, Cells = Cells, Divisions = Divisions))
  }
  

import.processedData <- function (File, Image1, Image2, leafShape, meanAngle)
  {
  # This function imports the processed data from PointTracker
  # computed between two png images 'Image1' and 'Image2',
  # as well as the raw data for all tracked images.
  # The xls file should be in a subfolder of the working directory called 'Growth'.
  #
  # Below is a description of the arguments:
  #
  #   - 'File' is the xls file containing raw and processed data. Use for example:
  #         File = "growth.xls"
  #
  #   - 'Image1' and 'Image2' are the png files for the first (initial)
  #     and second (final) image, respectively. Use for example:
  #         Image1 = "TL00.png"
  #         Image2 = "TL10.png"
  #
  #   - Set 'leafShape' to TRUE if the leaf outlines are available as txt files
  #     saved as XY coordinates using ImageJ under 'Leaf Shape\Name_Of_Image_i.txt'.
  #     It will be used here only to compute 'meanAngle' (see below).
  #
  #   - 'meanAngle' is the angle by which the leaf will be rotated so that the leaf
  #     will be oriented upward. It should be given in degrees from the y-axis
  #     (once the 'AngleShift' from PointTracker has been taken into account).
  #     See also 'get.meanAngle()'.
  #     It will be used only to compute Kpar and Kper (see below).


  # Image dimensions
  image_properties1 <- GDALinfo(paste("Processed", Image1, sep = "/"), silent = T)
  wd1 <- as.numeric(image_properties1["columns"]) # width of Image1 (pixels)
  ht1 <- as.numeric(image_properties1["rows"]) # height of Image1 (pixels)
  image_properties2 <- GDALinfo(paste("Processed", Image2, sep = "/"), silent = T)
  wd2 <- as.numeric(image_properties2["columns"]) # width of Image2 (pixels)
  ht2 <- as.numeric(image_properties2["rows"]) # height of Image2 (pixels)
  Image.dim <- data.frame(WidthImage1 = wd1, HeightImage1 = ht1, WidthImage2 = wd2, HeightImage2 = ht2)


  # Import growth file
  dataImport <- read.table(paste("Growth", File, sep = "/"), sep = ",", header = F, fill = T, blank.lines.skip = F, col.names = 1:1000)
  activeCol <- vector()
  for (Col in 1:ncol(dataImport)) { activeCol <- c(activeCol, length(which(is.na(dataImport[Col])))) }
  if (length(which(activeCol == nrow(dataImport))) != 0) { dataImport <- dataImport[-which(activeCol == nrow(dataImport))] }


  # Image properties stored in the growth file
  Images <- dataImport[which(dataImport[1] == "List of the images used"), 2:ncol(dataImport)]
  Images <- Images[!is.na(Images) & Images != ""]
  if (which(Image2 == Images) != 1 + which(Image1 == Images)) { stop("Image2 should be immediately posterior to Image1.", call. = F) }

  Info <- dataImport[which(dataImport[1] == "Data") + 2:5, 1:length(which(!is.na(dataImport[(which(dataImport[1] == "Data") + 5), ])))]
  InfoVertices <- data.frame(Images = NULL, ShiftX = NULL, ShiftY = NULL, AngleShift = NULL, Time = NULL, ScalingX = NULL, ScalingY = NULL)
  for (imag in seq(2, ncol(Info), 2)) { InfoVertices <- rbind(InfoVertices, data.frame(Images = as.character(Info[1, imag]), ShiftX = as.numeric(as.character(Info[2, imag])), ShiftY = as.numeric(as.character(Info[2, imag+1])), AngleShift = as.numeric(as.character(Info[3, imag])), Time = as.numeric(as.character(Info[3, imag+1])), ScalingX = as.numeric(as.character(Info[4, imag])), ScalingY = as.numeric(as.character(Info[4, imag+1])))) }
  InfoVertices$Images <- as.factor(InfoVertices$Images)

  if (missing(leafShape))
    {
    ifelse (file.exists(paste("Leaf Shape/", substr(Image1, 1, nchar(Image1) - 3), "txt", sep = "")) &
            file.exists(paste("Leaf Shape/", substr(Image2, 1, nchar(Image2) - 3), "txt", sep = "")),
            leafShape <- T, leafShape <- F)
    }
  if (missing(meanAngle))
    {
    ifelse (leafShape, meanAngle <- get.meanAngle(InfoVertices), meanAngle <- 0)
    }
  
  GrowthComputationParameters <- dataImport[which(dataImport[1] == "Estimation method"), 2:3]
  names(GrowthComputationParameters) <- c("EstimationMethod", "NbPoints")
  GrowthComputationParameters$EstimationMethod <- as.character(GrowthComputationParameters$EstimationMethod)
  GrowthComputationParameters$NbPoints <- as.integer(as.character(GrowthComputationParameters$NbPoints))

    
  # Growth
  Growth <- dataImport[(which(dataImport[1] == "Growth per image") + 2):(which(dataImport[1] == "Actual cell shapes") - 1), 1:10]
  is.Image.start <- which(Growth[1] != "")
  is.Image2.start <- which(Growth[1] == Image2)
  ifelse (is.Image2.start != is.Image.start[length(is.Image.start)], is.Image2.end <- is.Image.start[which(is.Image.start == is.Image2.start) + 1] - 1, is.Image2.end <- nrow(Growth))
  GrowthWall <- Growth[(is.Image2.start+1):is.Image2.end, -(1:8)]
  GrowthWall[, 1] <- factor(GrowthWall[, 1])
  GrowthWall[, 2] <- as.numeric(as.character(GrowthWall[, 2]))
  Growth <- Growth[(is.Image2.start+1):is.Image2.end, -c(1, 8:10)]
  for (Col in 1:ncol(Growth)) { names(Growth)[Col] <- as.character(dataImport[(which(dataImport[1] == "Growth per image") + 1), 1 + Col]) }
  if (length(which(Growth$Cell == "")) != 0) { Growth <- Growth[-which(Growth$Cell == ""), ] }
  Growth[, 1] <- factor(Growth[, 1])
  for (Col in 2:6) { Growth[, Col] <- as.numeric(as.character(Growth[, Col])) }
  Growth$Anisotropy <- 1 - Growth[, 4]/Growth[, 3]
  if (length(which(is.na(Growth[, 3]))) > 0) { Growth$Anisotropy[which(is.na(Growth[, 3]))] <- 0 }


  # Compute the growth rates parallel (ky, or kml) and perpendicular (kx) to the midvein
  Growth$kx <- NA
  Growth$ky <- NA
  vx <- c(cos(meanAngle * pi/180), sin(meanAngle * pi/180))               # the X unit vector (parallel to leaf width)    #ifelse (!HorizontalToVertical, vx <- c(1, 0), vx <- c(0, 1))
  vy <- c(cos(meanAngle * pi/180 + pi/2), sin(meanAngle * pi/180 + pi/2)) # the Y unit vector (parallel to midvein)       #ifelse (!HorizontalToVertical, vy <- c(0, 1), vy <- c(1, 0))
  for (i in 1:nrow(Growth))
    {
    GT <- paramsToTensor(list(kmaj = Growth[i, "kmaj (1/h)"],
                              kmin = Growth[i, "kmin (1/h)"],
                              theta = Growth[i, "Orientation of the major axis (Â° to the x axis)"] * pi/180,
                              phi = Growth[i, "phi (1/h)"]))
    Growth$kx[i] <- as.numeric(t(vx) %*% GT %*% vx) # growth rate parallel to the X axis (leaf width)
    Growth$ky[i] <- as.numeric(t(vy) %*% GT %*% vy) # growth rate parallel to the Y axis (midvein)
    }


  # Densified shapes
  Shapes <- dataImport[(which(dataImport[1] == "Actual cell shapes") + 2):(which(dataImport[1] == "Data") - 2), 1:ncol(dataImport)]
  is.Image.start <- which(Shapes[1] != "")
  is.Image2.start <- which(Shapes[1] == Image2)
  ifelse (is.Image2.start != is.Image.start[length(is.Image.start)], is.Image2.end <- is.Image.start[which(is.Image.start == is.Image2.start) + 1] - 1, is.Image2.end <- nrow(Shapes))
  Shapes <- Shapes[(is.Image2.start+1):is.Image2.end, -1]
  names(Shapes)[1] <- "Cell"
  names(Shapes)[2] <- "Begin/End"
  for (Col in seq(3, ncol(Shapes), 2)) { names(Shapes)[Col] <- paste("X", (Col-1)/2, sep = "") }
  for (Col in seq(4, ncol(Shapes), 2)) { names(Shapes)[Col] <- paste("Y", (Col-2)/2, sep = "") }
  Shapes[, 1] <- factor(Shapes[, 1])
  Shapes[, 2] <- factor(Shapes[, 2])
  for (Col in 3:ncol(Shapes)) { Shapes[, Col] <- as.numeric(as.character(Shapes[, Col])) }


  # Vertices coordinates
  Vertices <- dataImport[(which(dataImport[1] == "Data") + 6):(which(dataImport[1] == "Cells") - 1), 1:length(which(!is.na(dataImport[(which(dataImport[1] == "Data") + 6), ])))]
  names(Vertices)[1] <- "Point"
  for (Col in seq(2, ncol(Vertices), 2)) { names(Vertices)[c(Col, Col+1)] <- as.character(dataImport[(which(dataImport[1] == "Data") + 2), Col]) }
  Vertices[, 1] <- factor(Vertices[, 1])
  for (Col in 2:ncol(Vertices)) { Vertices[, Col] <- as.numeric(as.character(Vertices[, Col])) }


  # Cells to vertices
  Cells <- dataImport[(which(dataImport[1] == "Cells") + 1):(which(dataImport[1] == "Divisions") - 1), ]
  Cells[, 1] <- factor(Cells[, 1])
  for (Col in 2:ncol(Cells)) { Cells[, Col] <- as.numeric(as.character(Cells[, Col])) }
  activeCol <- vector()
  for (Col in 1:ncol(Cells)) { activeCol <- c(activeCol, length(which(is.na(Cells[Col])))) }
  if (length(which(activeCol == nrow(Cells))) != 0) { Cells <- Cells[-which(activeCol == nrow(Cells))] }
  names(Cells)[1] <- "Cell"
  names(Cells)[2:ncol(Cells)] <- paste("Vertex", 2:ncol(Cells) - 1, sep = "")


  # Divisions
  if (which(dataImport[1] == "Divisions") < nrow(dataImport))
    {
    Divisions <- dataImport[(which(dataImport[1] == "Divisions") + 1):nrow(dataImport), 1:6]
    names(Divisions)[1] <- "Cell"
    for (Col in 2:ncol(Divisions)) { Divisions[, Col] <- as.numeric(as.character(Divisions[, Col])) }
    for (Col in 2:ncol(Divisions)) { names(Divisions)[Col] <- as.character(dataImport[(which(dataImport[1] == "Divisions")), Col]) }
    for (t in as.numeric(as.character(levels(as.factor(Divisions$Time)))))
      {
      Divisions$TimeHours1[Divisions$Time == t] <- InfoVertices$Time[t]    # Time before cell division (hours)
      Divisions$TimeHours2[Divisions$Time == t] <- InfoVertices$Time[t+1]  # Time after cell division (hours)
      }
    Divisions$Cell <- factor(Divisions$Cell)
    }
  else
    {
    Divisions = NULL
    }


  # Compute cell areas at Image1 and Image2
  Growth$AreaAtTime1 <- NA
  Growth$AreaAtTime2 <- NA
  for (i in 1:nrow(Growth))
    {
    cell <- as.character(Growth$Cell[i])
    vert <- na.omit(as.numeric(Cells[Cells$Cell == cell, -1]))
    c(Data1, pts1) := cell.coord.for.alignCells(Image1, vert, Vertices)
    c(Data2, pts2) := cell.coord.for.alignCells(Image2, vert, Vertices)
    c(aligned_pts, aligned_new_pts, Ps, Qs) := alignCells(cell, pts1, pts2, Data1, Data2, GrowthComputationParameters$NbPoints, Image1, Image2)
    Growth$AreaAtTime1[i] <- areapl(Ps) * 1e12
    Growth$AreaAtTime2[i] <- areapl(Qs) * 1e12
    }
    
    
  # Final result to be returned
  return (list(Image.dim = Image.dim, InfoVertices = InfoVertices, meanAngle = meanAngle,
               Shapes = Shapes, Vertices = Vertices, Cells = Cells, Divisions = Divisions,
               GrowthComputationParameters = GrowthComputationParameters, Growth = Growth, GrowthWall = GrowthWall))
  }



################################################################################
#                  Compute growth, dense shapes, and divisions                 #
################################################################################

                           
recompute.GrowthAndShapes <- function (Cells, Divisions, Vertices, InfoVertices, Image1, Image2, meanAngle, method = "BackwardDense", nb_points = 100, exp_correction = T)
  {
  # This function recomputes growth and dense cell shapes from raw data.
  
  # Check method used to compute growth
  ifelse (method %in% c("Forward", "ForwardDense"), at_start <- T, ifelse (method %in% c("Backward", "BackwardDense"), at_start <- F, stop(paste("'", method, "' tracking method is not implemented for this script.", sep = ""))))
  if (!method %in% c("ForwardDense", "BackwardDense")) { stop(paste("'", method, "' tracking method is not implemented for this script.", sep = "")) }

  # Define x and y axes using the leaf as the reference (y is always parallel to the midvein)
  vx <- c(cos(meanAngle * pi/180), sin(meanAngle * pi/180))               # the X unit vector (parallel to leaf width)
  vy <- c(cos(meanAngle * pi/180 + pi/2), sin(meanAngle * pi/180 + pi/2)) # the Y unit vector (parallel to midvein)

  # Prepare dataframes
  Growth <- data.frame(V1 = NULL, V2 = NULL, V3 = NULL, V4 = NULL, V5 = NULL, V6 = NULL, V7 = NULL, V8 = NULL, V9 = NULL, V10 = NULL, V11 = NULL)
  Shapes <- data.frame(Cell = NULL, V2 = NULL)
  
  # Image time
  time1 <- InfoVertices$Time[InfoVertices$Images == Image1]
  time2 <- InfoVertices$Time[InfoVertices$Images == Image2]
  DT <- time2 - time1
  
  # Find the current cells at the image considered
  levels.Cells <- find.Cells(Image1, Cells, Divisions, Vertices, InfoVertices)
  
  # Loop on each cell
  for (cell in levels.Cells)
    {
    # Find vertices
    vert <- na.omit(as.numeric(Cells[Cells$Cell == cell, -1]))
    
    # Compute dense shapes
    c(Data1, pts1) := cell.coord.for.alignCells(Image1, vert, Vertices)
    c(Data2, pts2) := cell.coord.for.alignCells(Image2, vert, Vertices)
    c(aligned_pts, aligned_new_pts, Ps, Qs) := alignCells(cell, pts1, pts2, Data1, Data2, nb_points, Image1, Image2)
        
    # Compute and save growth
    if (Image1 != Image2)
      {
      GT <- suppressWarnings(growthTensor(P = Ps, Q = Qs, DT = DT, at_start = at_start, exp_correction = exp_correction))
      if (is.nan(GT[1])) { stop(paste("\n\n", paste(rep("#", nchar(paste("# check", cell, "in", Image1, "and", Image2, "#"))), sep = "", collapse = ""),
                                      "\n # check", cell, "in", Image1, "and", Image2, "#\n",
                                      paste(rep("#", nchar(paste("# check", cell, "in", Image1, "and", Image2, "#"))), sep = "", collapse = ""), "\n"), call. = F) }
      c(kmaj, kmin, theta, phi) := tensorToParams(GT)
      kx <- as.numeric(t(vx) %*% GT %*% vx) # growth rate parallel to the X axis (leaf width)
      ky <- as.numeric(t(vy) %*% GT %*% vy) # growth rate parallel to the Y axis (midvein)
      #karea <- kmaj + kmin # = kx + ky
      a1 <- areapl(Ps)
      a2 <- areapl(Qs)
      karea <- log(a2/a1)/DT # This is the way 'karea' is implemented in PointTracker
      ifelse (kmaj == 0, aniso <- 0, aniso <- 1 - kmin/kmaj)
      }
    else
      {
      kmaj <- kmin <- theta <- phi <- kx <- ky <- a1 <- a2 <- karea <- aniso <- NA
      }
    Growth <- rbind(Growth, data.frame(V1 = cell, V2 = karea, V3 = kmaj, V4 = kmin, V5 = theta * 180/pi, V6 = phi, V7 = aniso, V8 = kx, V9 = ky, V10 = a1 * 1e12, V11 = a2 * 1e12))
    
    # Save dense shapes
    actual.nb_points <- nrow(Ps) # close to 'nb_points', but generally not equal to
    if (nrow(Shapes) > 0)
      {
      while (actual.nb_points > (ncol(Shapes)-2)/2)
        {
        Shapes <- cbind(Shapes, NA)
        names(Shapes)[ncol(Shapes)] <- paste("X", ncol(Shapes)-2, sep = "")
        }
      while (actual.nb_points < (ncol(Shapes)-2)/2)
        {
        Ps <- rbind(Ps, NA)
        Qs <- rbind(Qs, NA)
        actual.nb_points <- nrow(Ps)
        }
      }
    Shapes <- rbind(Shapes, data.frame(Cell = cell, V2 = "Begin", t(as.numeric(Ps)[order(rep(1:actual.nb_points, 2))])))
    Shapes <- rbind(Shapes, data.frame(Cell = cell, V2 = "End", t(as.numeric(Qs)[order(rep(1:actual.nb_points, 2))])))   
    }  

  # Finalize data frames
  names(Growth) <- c("Cell", "karea (1/h)", "kmaj (1/h)", "kmin (1/h)", "Orientation of the major axis (Â° to the x axis)", "phi (1/h)", "Anisotropy", "kx", "ky", "AreaAtTime1", "AreaAtTime2")
  #Growth$Cell <- as.factor(Growth$Cell)
  names(Shapes)[2] <- "Begin/End"
  names(Shapes)[seq(3, ncol(Shapes), 2)] <- paste("X", 1:((ncol(Shapes)-2)/2), sep = "")
  names(Shapes)[seq(4, ncol(Shapes), 2)] <- paste("Y", 1:((ncol(Shapes)-2)/2), sep = "")  
  #Shapes[, 1:2] <- as.factor(Shapes[, 1:2])
  
  return (list(Growth = Growth, Shapes = Shapes))
  }


process.AllDiv <- function (InfoVertices, Vertices, Cells, Divisions, meanAngle,
                            alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary,
                            optimal.interval = 24, leafShape, CSV = T, File, ini = 1)
  {
  # For each cell in 'Divisions$Cells', this function appends to the dataframe 'Divisions'
  # with the following information, when available:
  #
  #   - ID of the direct mother, birth time and cell cycle duration
  #
  #   - area at division, area of the daughters and symmetry of the division
  #
  #   - growth parameters computed for each cell during a period 'DT_Used_To_Compute_Growth',
  #     which ends at the time of actual division 'Divisions$TimeHours2'
  #     and starts at about 'optimal.interval' hours before the actual division.
  #     If 'Divisions$TimeHours2 - optimal.interval' is anterior to the first tracked image,
  #     then growth will be computed between (1) the first image of the dataset and (2) the image
  #     posterior to 'Divisions$TimeHours2' which is at closest to 'optimal.interval'.
  #
  # The dataframe will be saved as 'Divisions__<name of PointTracker csv file>.csv' in a subfolder 'Divisions',
  # only if 'CSV' is TRUE (the default).
  #
  # The argument 'ini' is the index of the image from which divisions start to be processed.
  # This is useful to deal with datasets in which some cells are not visible from the beginning;
  # however, it will introduce bias in the computation of growth in all cells that divide from the image 'ini' to 'optimal.interval',
  # even those that were visible from the beginning (because it is not possible to determine which are the new cells, since the ones
  # that share more than two points with a pre-existing cell will have a non-null area before 'ini').

  # --> edit 2016/09/22 #
  InfoVertices$Time <- round(InfoVertices$Time, 4)
  Divisions$TimeHours1 <- round(Divisions$TimeHours1, 4)
  Divisions$TimeHours2 <- round(Divisions$TimeHours2, 4)
  # edit 2016/09/22 <-- #

  # Prepare the data frame
  Divisions$Lineage <- NA
  Divisions$DirectMother <- NA
  Divisions$OldestAncestor <- NA
  Divisions$TimeHours1_Birth <- NA
  Divisions$TimeHours2_Birth <- NA
  Divisions$Area_At_Division <- NA
  Divisions$Area_Daughter1_At_Division <- NA
  Divisions$Area_Daughter2_At_Division <- NA
  Divisions$Symmetry <- NA
  Divisions$DT_Used_To_Compute_Growth <- NA
  Divisions$karea <- NA
  Divisions$kmaj <- NA
  Divisions$kmin <- NA
  Divisions$theta <- NA
  Divisions$phi <- NA
  Divisions$aniso <- NA
  Divisions$kx <- NA
  Divisions$ky <- NA
  
  # Cells at first image (to get the lineage number)
  cellsAtFirstImage <- find.Cells(InfoVertices$Images[ini], Cells, Divisions, Vertices, InfoVertices)
    
  # To save time, loop on each actual time of division
  tini <- InfoVertices$Time[ini]
  for (t2 in unique(Divisions$TimeHours2)[unique(Divisions$TimeHours2) > tini])
    {
    idx <- which(Divisions$TimeHours2 == t2)
    
    # Images at division
    t1 <- unique(Divisions$TimeHours1[idx])
    Image1_Division <- as.character(InfoVertices$Image[InfoVertices$Time == t1])
    Image2_Division <- as.character(InfoVertices$Image[InfoVertices$Time == t2])

    # Import leaf outline, cell shapes and vertices at 'Image1_Division' and 'Image2_Division'
    c(Growth_NotUsed, Shapes) := recompute.GrowthAndShapes(Cells, Divisions, Vertices, InfoVertices, Image1_Division, Image2_Division, meanAngle)
    c(Vert1X, Vert1Y, Shapes1X, Shapes1Y, LeafShape1, LeafSpline1) := scale.Objects(Image1_Division, Image2_Division, before = T, meanAngle, leafShape, alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary, Shapes, Cells, Vertices, InfoVertices)
    c(Vert2X, Vert2Y, Shapes2X, Shapes2Y, LeafShape2, LeafSpline2) := scale.Objects(Image1_Division, Image2_Division, before = F, meanAngle, leafShape, alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary, Shapes, Cells, Vertices, InfoVertices)

    # Loop on each division event at t2
    for (i in idx)
      {
      # Cell ID
      cell <- as.character(Divisions$Cell[i])
      
      # Lineage number
      c(cellAtFirstImage, direct.mother) := find.OldestAncestor(cell, InfoVertices$Time[ini], Divisions)
      Lineage <- which(cellAtFirstImage == cellsAtFirstImage)
      
      # Birth time
      if (length(direct.mother) == 0)
        {
        TimeHours2_Birth <- TimeHours1_Birth <- direct.mother <- NA
        }
      else
        {
        Image1_Birth <- as.character(InfoVertices$Image[InfoVertices$Time == Divisions$TimeHours1[Divisions$Cell == direct.mother]])
        TimeHours1_Birth <- InfoVertices$Time[InfoVertices$Image == Image1_Birth]
        Image2_Birth <- as.character(InfoVertices$Image[InfoVertices$Time == Divisions$TimeHours2[Divisions$Cell == direct.mother]])
        TimeHours2_Birth <- InfoVertices$Time[InfoVertices$Image == Image2_Birth]
        }
  
      # Daughter cells ID
      Daughter1 <- paste("Cell", Divisions[i, "Daughter cell 1"])
      Daughter2 <- paste("Cell", Divisions[i, "Daughter cell 2"])
  
      # ID of the vertices forming the mother and daughter cells
      vert_mother <- na.omit(as.numeric(Cells[Cells$Cell == cell, -1]))
      vert_D1 <- na.omit(as.numeric(Cells[Cells$Cell == Daughter1, -1]))
      vert_D2 <- na.omit(as.numeric(Cells[Cells$Cell == Daughter2, -1]))

      # Coordinates of the vertices of the mother and daughter cells at division
      XY1_mother <- cell.coord(vert_mother, Vert1X, Vert1Y, Vertices) # at 'Image1_Division'
      XY2_mother <- cell.coord(vert_mother, Vert2X, Vert2Y, Vertices) # at 'Image2_Division'
      XY2_D1 <- cell.coord(vert_D1, Vert2X, Vert2Y, Vertices)
      XY2_D2 <- cell.coord(vert_D2, Vert2X, Vert2Y, Vertices)
    
      # Area of the cells
      Area_At_Division <- mean(areapl(na.omit(XY1_mother)), areapl(na.omit(XY2_mother)))
      Area_Daughters <- areapl(na.omit(XY2_D1)) + areapl(na.omit(XY2_D2)) # = areapl(na.omit(XY2_mother))
      Area_Daughter1_At_Division <- areapl(na.omit(XY2_D1)) * Area_At_Division / Area_Daughters  # to correct for the growth between 'Image1_Division' and 'Image2_Division'
      Area_Daughter2_At_Division <- areapl(na.omit(XY2_D2)) * Area_At_Division / Area_Daughters  # to correct for the growth between 'Image1_Division' and 'Image2_Division'
  
      # Search images and oldest ancestor to compute growth
    
        # 'Image2' is 'Image2_Division'
        Image2 <- Image2_Division
        Time2 <- t2
        
        # 'Image1' is image whose time minimises the time of 'Image2' minus 'optimal.interval'
        ID_Closest_Image_To_Optimal_Interval <- which.min(((Time2 - optimal.interval) - InfoVertices$Time[ini:nrow(InfoVertices)])^2) + (ini - 1)
        Image1 <- as.character(InfoVertices$Images[ID_Closest_Image_To_Optimal_Interval])
        Time1 <- InfoVertices$Time[ID_Closest_Image_To_Optimal_Interval]
        c(oldest.ancestor, dm) := find.OldestAncestor(cell, Time1, Divisions)
        
        # If the beginning of the dataset has been reached by 'Image1', extend the period after 'Image2' to compute growth
        if (Image1 == InfoVertices$Images[ini])
          {
          ID_Closest_Image_To_Optimal_Interval <- which.min(((Time1 + optimal.interval) - InfoVertices$Time[ini:nrow(InfoVertices)])^2) + (ini - 1)
          Image2 <- as.character(InfoVertices$Images[ID_Closest_Image_To_Optimal_Interval])
          Time2 <- InfoVertices$Time[ID_Closest_Image_To_Optimal_Interval]
          }
  
      # Compute growth
      c(Growth, Shapes_NotUsed) := recompute.GrowthAndShapes(Cells, Divisions, Vertices, InfoVertices, Image1, Image2, meanAngle)
      c(karea, kmaj, kmin, theta, phi, aniso, kx, ky) := Growth[Growth$Cell == oldest.ancestor, 2:9]
      
      # Assign in the dataframe
      Divisions$Lineage[i] <- Lineage
      Divisions$DirectMother[i] <- direct.mother
      Divisions$OldestAncestor[i] <- cellAtFirstImage
      Divisions$TimeHours1_Birth[i] <- TimeHours1_Birth
      Divisions$TimeHours2_Birth[i] <- TimeHours2_Birth
      Divisions$Area_At_Division[i] <- Area_At_Division
      Divisions$Area_Daughter1_At_Division[i] <- Area_Daughter1_At_Division
      Divisions$Area_Daughter2_At_Division[i] <- Area_Daughter2_At_Division
      Divisions$Symmetry[i] <- min(c(Area_Daughter1_At_Division, Area_Daughter2_At_Division)) / (Area_At_Division/2)
      Divisions$DT_Used_To_Compute_Growth[i] <- Time2 - Time1 # should be close to 'optimal.interval'
      Divisions$karea[i] <- karea
      Divisions$kmaj[i] <- kmaj
      Divisions$kmin[i] <- kmin
      Divisions$theta[i] <- theta
      Divisions$phi[i] <- phi
      Divisions$aniso[i] <- aniso
      Divisions$kx[i] <- kx
      Divisions$ky[i] <- ky
      }
    }

  # Compute cell cycle duration and return 'Divisions'
  Divisions$Mean_Time_Division <- (Divisions$TimeHours1 + Divisions$TimeHours2)/2
  Divisions$Mean_Time_Birth <- (Divisions$TimeHours1_Birth + Divisions$TimeHours2_Birth)/2
  Divisions$Cell_Cycle_Duration <- Divisions$Mean_Time_Division - Divisions$Mean_Time_Birth
  
  # Save in csv file and return
  if (CSV)
    {
    if (!file.exists("Divisions")) { dir.create("Divisions") } 
    write.csv(Divisions, paste("Divisions/Divisions__", File, sep = ""), row.names = F)
    }
  return(Divisions = Divisions)
  }

                        
process.DivWithInterval <- function (InfoVertices, Vertices, Cells, Divisions, meanAngle, Growth, Shapes,
                                     alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary,
                                     Image1, Image2, leafShape, before, File, ini)
  {
  # This function returns a dataframe 'Div' in which each cell at the image considered
  # ('Image1' if before is TRUE, 'Image2' if before is FALSE) is associated with
  # its direct mother (if any), its oldest ancestor at 'Image1' (if any),
  # its status of division, its status of competence for division,
  # its division class as defined by the user (not implemented yet) and its growth parameters,
  # all computed between 'Image1' and 'Image2'. Division parameters (area at division,
  # symmetry, cell cycle duration if date of birth available), as well as the blade dimensions
  # and cell areas at 'Image1' and 'Image2', are also returned.
  #
  # It makes more sense to set 'before' to FALSE, although TRUE will also work
  # (but if some cells have undergone several rounds of division, it will not appear in 'Div').

  # --> edit 2016/09/22 #
  InfoVertices$Time <- round(InfoVertices$Time, 4)
  # edit 2016/09/22 <-- #

  # Import leaf outline, cell shapes and vertices
  if (missing(leafShape))
    {
    ifelse (file.exists(paste("Leaf Shape/", substr(Image1, 1, nchar(Image1) - 3), "txt", sep = "")) &
            file.exists(paste("Leaf Shape/", substr(Image2, 1, nchar(Image2) - 3), "txt", sep = "")),
            leafShape <- T, leafShape <- F)
    }
  c(Vert1X, Vert1Y, Shapes1X, Shapes1Y, LeafShape1, LeafSpline1) := scale.Objects(Image1, Image2, before = T, meanAngle, leafShape, alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary, Shapes, Cells, Vertices, InfoVertices)
  c(Vert2X, Vert2Y, Shapes2X, Shapes2Y, LeafShape2, LeafSpline2) := scale.Objects(Image1, Image2, before = F, meanAngle, leafShape, alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary, Shapes, Cells, Vertices, InfoVertices)

  # Retrieve processed divisions or process all divisions
  if (file.exists(paste("Divisions/Divisions__", File, sep = "")))
    {
    Divisions <- read.csv(paste("Divisions/Divisions__", File, sep = ""))
    names(Divisions) <- gsub("\\.", " ", names(Divisions))
    # --> edit 2016/09/22 #
    Divisions$TimeHours1 <- round(Divisions$TimeHours1, 4)
    Divisions$TimeHours2 <- round(Divisions$TimeHours2, 4)
    # edit 2016/09/22 <-- #
    }
  else if (!is.null(Divisions))
    {
    Divisions <- process.AllDiv(InfoVertices, Vertices, Cells, Divisions, meanAngle,
                                alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary,
                                optimal.interval = 24, leafShape, CSV = T, File, ini)
    # --> edit 2016/09/22 #
    Divisions$TimeHours1 <- round(Divisions$TimeHours1, 4)
    Divisions$TimeHours2 <- round(Divisions$TimeHours2, 4)
    # edit 2016/09/22 <-- #
    }
      
  # Blade dimensions  
  if (leafShape)
    {
    bladeLength1 <- max(LeafSpline1$y)
    bladeWidth1 <- max(LeafSpline1$x) - min(LeafSpline1$x)
    bladeLength2 <- max(LeafSpline2$y)
    bladeWidth2 <- max(LeafSpline2$x) - min(LeafSpline2$x)
    }
  else
    {
    bladeWidth2 <- bladeWidth1 <- bladeLength2 <- bladeLength1 <- NA
    }
  
  # Find the current cells at the image considered
  ifelse (before, Image <- Image1, Image <- Image2)
  levels.Cells <- find.Cells(Image, Cells, Divisions, Vertices, InfoVertices)
  
  # Search divided cells and competent cells.
  # Divided cells are those who have divided between time1 and time2.
  # Competent cells are those who will divide again after time2;
  # a cell who has just divided is NOT considered as competent at time2,
  # except if one of its daughters divides again after time2.
  time1 <- InfoVertices$Time[InfoVertices$Images == Image1]
  time2 <- InfoVertices$Time[InfoVertices$Images == Image2]
  has.Divided <- Divisions[Divisions$TimeHours2 > time1 & Divisions$TimeHours2 <= time2, ] # include cells which have undergone several rounds of divisions
  mothers <- as.character(has.Divided$Cell)
  daughters <- paste("Cell", c(has.Divided[, "Daughter cell 1"], has.Divided[, "Daughter cell 2"]))
  competent.cells <- as.character(Divisions$Cell[Divisions$TimeHours2 > time2 & Divisions$Cell %in% levels.Cells])
  if (before)
    {
    divided.cells <- mothers
    for (m in 1:length(mothers))
      {
      competent.daughters <- which(Divisions$TimeHours2 > time2 & Divisions$Cell %in% daughters[c(m, length(daughters)/2+m)])
      if (length(competent.daughters) > 0) # add mothers with competent daughters to the competent cells
        {
        m.AtImage1 <- mothers[m]
        while (!m.AtImage1 %in% levels.Cells)
          {
          m.AtImage1 <- rep(mothers, 2)[which(daughters == m.AtImage1)] # find ancestor at Image1
          }
        competent.cells <- c(competent.cells, m.AtImage1)
        }
      }
    competent.cells <- unique(competent.cells) # because ancestors may have been incorporated several times
    }
  else
    {
    divided.cells <- daughters
    }

  # Load the user table of division classes (if any)
  DivType <- NULL
  if (file.exists("Divisions/DivisionClass.xlsx"))
    {
    if (!"XLConnect" %in% installed.packages()[, "Package"]) { install.packages("XLConnect") }
    if (!"package:XLConnect" %in% search()) { library("XLConnect", character.only = T) } # the package 'XLConnect' cannot be loaded at the beginning as it impairs 'rJava'
    workBook <- loadWorkbook("Divisions/DivisionClass.xlsx")
    sheetNames <- getSheets(workBook)
    if ("DivClass" %in% sheetNames)
      {
      DivType <- readWorksheet(workBook, sheet = "DivClass")
      DivType$Cell <- as.factor(DivType$Cell)
      if (nrow(DivType) != nrow(Divisions)) { stop("The number of divisions in the sheet 'DivClass' do not match the number of divisions in PointTracker") }
      }
    }

  # Prepare the dataframe
  Div <- data.frame(Image1 = NULL, Image2 = NULL, Time1 = NULL, Time2 = NULL,
                    BladeLength1 = NULL, BladeLength2 = NULL, BladeWidth1 = NULL, BladeWidth2 =NULL,
                    Cell = NULL, DirectMother = NULL, OldestAncestorAtImage1 = NULL,
                    HasDivided = NULL, IsCompetent = NULL, DivisionClass = NULL, CompetenceClass = NULL,
                    CompetenceClassDaughter1 = NULL, CompetenceClassDaughter2 = NULL, # edit 2016/09/22
                    Xm1 = NULL, Ym1 = NULL, Xm2 = NULL, Ym2 = NULL,
                    karea = NULL, kmaj = NULL, kmin = NULL, theta = NULL, phi = NULL, aniso = NULL, kx = NULL, ky = NULL,
                    AreaAtTime1 = NULL, AreaAtTime2 = NULL,
                    Area_At_Division = NULL, Symmetry = NULL, Cell_Cycle_Duration = NULL, Number_Of_Divisions = NULL) # edit 2015/12/15
          
  # Loop on each cell
  for (cell in levels.Cells)
    {
    vert <- na.omit(as.numeric(Cells[Cells$Cell == cell, -1]))
    classD2 <- classD1 <- comp.class <- div.class <- NA # edit 2016/06/22
    is.comp <- has.div <- F

    # Search for division and competence in the user table 'DivType'.
    #
    # Warning: which class will be saved if several rounds of division occured in the interval?
    # Because the loop is on 'levels.Cells':
    #
    #   - if 'before' is FALSE:
    #       - the division class will be the one of the LAST division;
    #       - the competence class will be the one of the division immediately posterior.
    #
    #   - if 'before' is TRUE:
    #       - the division class will be the one of the FIRST division;
    #       - the competence class will be the one of the division immediately posterior.
    #         If both daughter cells are competent and their competence class is not the same,
    #         THE COMPETENCE CLASS WITH THE SMALLEST INDEX WILL BE RETURNED.
    #         For example, if 'daughter 1' is competent for 'class 3' and
    #         'daughter 2' is competent for 'class 1', 'class 1' will be returned.

      # Search for division class 
      if (cell %in% divided.cells)
        {
        has.div <- T
        if (is.null(DivType)) { div.class <- 1 }
        else
          {
          ifelse (before, direct.mother <- cell, direct.mother <- as.character(rep(Divisions$Cell, 2)[which(paste("Cell", c(Divisions[, "Daughter cell 1"], Divisions[, "Daughter cell 2"])) == cell)]))
          div.class <- DivType$Class[DivType$Cell == direct.mother]
          }
        }
  
      # Search for competence class
      if (cell %in% competent.cells)
        {
        is.comp <- T
        if (is.null(DivType))
          {
          comp.class <- 1
          # --> edit 2016/09/22 #
          if (before & cell %in% divided.cells)
            {
            D1 <- paste("Cell", Divisions[Divisions$Cell == cell, "Daughter cell 1"])
            D2 <- paste("Cell", Divisions[Divisions$Cell == cell, "Daughter cell 2"])
            D1.1 <- Divisions[Divisions$Cell == D1, "Daughter cell 1"]
            if (length(D1.1) != 0) { classD1 <- 1 }
            D2.1 <- Divisions[Divisions$Cell == D2, "Daughter cell 1"]
            if (length(D2.1) != 0) { classD2 <- 1 }
            }
          # edit 2016/09/22 <-- #
          }
        else
          {
          if (before & cell %in% divided.cells)
            {
            D1 <- paste("Cell", Divisions[Divisions$Cell == cell, "Daughter cell 1"])
            D2 <- paste("Cell", Divisions[Divisions$Cell == cell, "Daughter cell 2"])
            classD1 <- DivType$Class[DivType$Cell == D1]
            ifelse (length(classD1) != 0, comp.class <- classD1, classD1 <- NA) # edit 2016/09/22
            classD2 <- DivType$Class[DivType$Cell == D2]
            if (length(classD2) != 0)
              {
              if (!is.na(classD1)) { comp.class <- min(c(classD1, classD2)) } # edit 2016/09/22
              else { comp.class <- classD2 }
              }
            else { classD2 <- NA } # edit 2016/09/22
            }
          else
            {
            comp.class <- DivType$Class[DivType$Cell == cell]
            }
          }
        }

    # --> edit 2016/09/22 #
    # If 'before' is FALSE and the cell has divided,
    # the oldest ancestor present at 'Image1' is used
    # to compute growth between 'Image1' and 'Image2'.
    c(oldest.ancestor, direct.mother) := find.OldestAncestor(cell, time1, Divisions)  # NB: oldest.ancestor == cell if before is TRUE or if the cell has not divided
    if (length(direct.mother) == 0) { direct.mother <- NA }
    if (!before & cell %in% divided.cells)
      {
      vert <- na.omit(as.numeric(Cells[Cells$Cell == oldest.ancestor, -1]))
      }
    # edit 2016/09/22 <-- #

    # Find centroids of the cell (or of its ancestor if the cell has divided).
    # Note that:
    #
    #   - Division status should be related to the MEAN centroid coordinates
    #     in 'Image1' and 'Image2'
    #
    #   - Competence status should be related to the centroid coordinates in
    #       - 'Image2' if 'before' is FALSE
    #       - 'Image1' if 'before' is TRUE
    #
    XY1 <- cell.coord(vert, Vert1X, Vert1Y, Vertices)
    c(Xm1, Ym1) := get.Centroid(XY1) 
    XY2 <- cell.coord(vert, Vert2X, Vert2Y, Vertices)
    c(Xm2, Ym2) := get.Centroid(XY2)
    
    # Retrieve growth parameters
    c(karea, kmaj, kmin, theta, phi, aniso, kx, ky, a1, a2) := Growth[Growth$Cell == oldest.ancestor, 2:11]
    
    # Retrieve division parameters
    if (cell %in% divided.cells)
      {
      if (before)
        {
        Area_At_Division <- Divisions$Area_At_Division[Divisions$Cell == cell]
        Symmetry <- Divisions$Symmetry[Divisions$Cell == cell]
        Cell_Cycle_Duration <- Divisions$Cell_Cycle_Duration[Divisions$Cell == cell]
        }
      else
        {
        Area_At_Division <- Divisions$Area_At_Division[Divisions$Cell == direct.mother]
        Symmetry <- Divisions$Symmetry[Divisions$Cell == direct.mother]
        Cell_Cycle_Duration <- Divisions$Cell_Cycle_Duration[Divisions$Cell == direct.mother]
        }
      n.div <- find.MaxNumberOfDivisions(cell, time1, time2, Divisions, before) # edit 2015/12/15
      }
    else
      {
      Area_At_Division <- Symmetry <- Cell_Cycle_Duration <- n.div <- NA # edit 2015/12/15
      }
      
    # Assignment in the dataframe
    Div <- rbind(Div, data.frame(Image1 = Image1, Image2 = Image2, Time1 = time1, Time2 = time2,
                                 BladeLength1 = bladeLength1, BladeLength2 = bladeLength2, BladeWidth1 = bladeWidth1, BladeWidth2 = bladeWidth2,
                                 Cell = cell, DirectMother = direct.mother, OldestAncestorAtImage1 = oldest.ancestor,
                                 HasDivided = has.div, IsCompetent = is.comp, DivisionClass = div.class, CompetenceClass = comp.class,
                                 CompetenceClassDaughter1 = classD1, CompetenceClassDaughter2 = classD2, # edit 2016/09/22
                                 Xm1 = Xm1, Ym1 = Ym1, Xm2 = Xm2, Ym2 = Ym2,
                                 karea = karea, kmaj = kmaj, kmin = kmin, theta = theta, phi = phi, aniso = aniso, kx = kx, ky = ky,
                                 AreaAtTime1 = a1, AreaAtTime2 = a2,
                                 Area_At_Division = Area_At_Division, Symmetry = Symmetry, Cell_Cycle_Duration = Cell_Cycle_Duration, Number_Of_Divisions = n.div)) # edit 2015/12/15
    }
  
  # Return 'Div' (NB: 'DivType' could be returned here as well)
  return (Div = Div)
  }



################################################################################
#                                   Plot maps                                  #
################################################################################


# Perform initial checks on parameters given by the user

initialChecks <- function (scaleBar, xlim, ylim, zlim, aniso.threshold, aniso.lwd,
                           n.pred, polynomDegree, rangeParameter, contour.levels,
                           leafShapeKriging, cellularScale, plotResidual, exportCellValues,
                           div.parameter, check.type)
  {
  # This function checks if parameters for the plotting functions are correct.
  
  if (check.type %in% c("growth", "kriging", "division"))
    {
    if (!missing(scaleBar))
      {
      if (length(scaleBar) != 1 | !is.numeric(scaleBar))
        {
        stop("The parameter 'scaleBar' should be a numeric values in microns.")
        }
      }
    
    if (!missing(xlim))
      {
      if (length(xlim) != 2 | !is.numeric(xlim) | length(which(is.na(xlim))) != 0)
        {
        stop("The parameter 'xlim' should be either missing or given as two numeric values.")
        }
      }
    
    if (!missing(ylim))
      {
      if (length(ylim) != 2 | !is.numeric(ylim) | length(which(is.na(ylim))) != 0)
        {
        stop("The parameter 'ylim' should be either missing or given as two numeric values.")
        }
      }
      
    if (!missing(zlim))
      {
      if (length(zlim) != 2 | !is.numeric(zlim) | length(which(is.na(zlim))) != 0)
        {
        stop("The parameter 'zlim' should be either missing or given as two numeric values.")
        }
      }
          
    # Specific to growth and surface estimation
    if (check.type %in% c("growth", "kriging"))
      {
      if (aniso.threshold < 0 | aniso.threshold > 1 | length(aniso.threshold) != 1)
        {
        stop("The parameter 'aniso.threshold' should lie between 0 and 1.")
        }
      
      if (!missing(aniso.lwd))
        {
        if (aniso.lwd <= 0 | length(aniso.lwd) != 1)
          {
          stop("The parameter 'aniso.lwd' should be a strictly positive value.")
          }
        }
      
      # Specific to kriging
      if (check.type == "kriging")
        {
        if (!missing(n.pred))
          {
          if (n.pred < 10 | n.pred > 1000 | length(n.pred) != 1 | as.integer(n.pred) != n.pred)
            {
            stop("The parameter 'n.pred' should be an integer between 10 and 1000.")
            }
          }
                   
        if (!missing(contour.levels))
          {  
          if (length(contour.levels) < 1 | !is.numeric(contour.levels) | length(which(is.na(contour.levels))) != 0)
            {
            stop("The parameter 'contour.levels' should be a numeric vector with at least 1 value.")
            }
          }
        
        if (!missing(polynomDegree))
          {
          if (polynomDegree < 2 | polynomDegree > 6 | length(polynomDegree) != 1 | as.integer(polynomDegree) != polynomDegree)
            {
            stop("The parameter 'polynomDegree' should be an integer between 2 and 6.")
            }    
          }
        if (!missing(rangeParameter))
          {
          if (rangeParameter <= 0 | length(rangeParameter) != 1)
            {
            stop("The parameter 'rangeParameter' should be a strictly positive number.")
            }
          }
        if (leafShapeKriging)
          {
          if (cellularScale | plotResidual | exportCellValues)
            {
            stop("Cell values are not compatible with kriging extrapolation over the whole leaf surface:
                  'cellularScale', 'plotResidual' and 'exportCellValues' should be FALSE when 'leafShapeKriging' is TRUE.")
            }
          }
        if (plotResidual)
          {
          if (!cellularScale)
            {
            stop("Plotting residuals is possible at the cellular level only:
                  'cellularScale' should be TRUE when 'plotResidual' is TRUE.")
            }
          }
        if (exportCellValues)
          {
          if (!cellularScale & !plotResidual)
            {
            stop("Exporting cell values is only possible if they are computed:
                  'cellularScale' or 'plotResidual' should be TRUE when 'exportCellValues' is TRUE.")
            }
          }
        }
      }

    # Specific to division
    else
      {
      #if (!div.parameter %in% c("div&comp" ,"div", "comp", "CellArea", "AreaAtTime1", "AreaAtTime2", "Area_At_Division", "Symmetry", "Cell_Cycle_Duration"))
      #  {
      #  stop("The argument 'div.parameter' should be one string within:\n\"div&comp\", \"div\", \"comp\", \"CellArea\", \"AreaAtTime1\", \"AreaAtTime2\", \"Area_At_Division\", \"Symmetry\" or \"Cell_Cycle_Duration\"")
      #  }
      # --> edit 2014/09/11 and 2015/12/15 #
      if (!div.parameter %in% c("div&comp" ,"div", "comp", "progeny", "CellArea", "AreaAtTime1", "AreaAtTime2", "Area_At_Division", "Symmetry", "Cell_Cycle_Duration", "Number_Of_Divisions"))
        {
        stop("The argument 'div.parameter' should be one string within:\n\"div&comp\", \"div\", \"comp\",  \"progeny\", \"CellArea\", \"AreaAtTime1\", \"AreaAtTime2\", \"Area_At_Division\", \"Symmetry\", \"Cell_Cycle_Duration\" or \"Number_Of_Divisions\"")
        }
      # edit 2014/09/11 and 2015/12/15 <-- #
      }
    }
    
  else
    {
    stop("The argument 'check.type' should be one string within \"growth\", \"kriging\" or \"division\"")
    }
  }


set.zlim <- function (k, Percent, round.zlim, zlim, fix.min, fix.max, colorPaletteOfGFtbox, n.colors = 1000, round.value = 0)
  {
  # This function sets the minimal and maximal growth values to be plotted
  # as a function of user's choices.
  #
  #   - 'k' is the input vector containing the growth values. Note that the function also returns
  #     a 'k' output vector, which takes 'fix.min' and 'fix.max' into account.
  #
  #   - Setting 'Percent' to TRUE will make growth to be plotted in %·h^-1 rather than in h^-1.
  #
  #   - Setting 'round.zlim' to TRUE (default) will make:
  #       - the maximum to be rounded to the multiple of 1%·h^-1 (or 0.01 h^-1)
  #         immediately superior to the maximum observed
  #       - the minimun to be 0, or the negative multiple of 1%·h^-1 (or 0.01 h^-1)
  #         immediately inferior to the minimum observed if this value is negative
  #
  #   - Alternatively, zlim can be given either in %·h^-1 or in h^-1, e.g. as, respectively:
  #       - zlim = c(0.5, 4.2)
  #       - zlim = c(0.005, 0.042)
  #
  #   - Setting 'fix.min' and/or 'fix.max' to TRUE will make the values below the minimum
  #     or above the maximum to be plotted with the color of the minimum or the maximum, respectively.
  #     Otherwise, the cell will be plotted with no color (only the cell outline will be visible).
  # 
  # This function also returns:
  #
  #   - the values of ticks to be plotted on the growth scale, 'ticksGrowthScale',
  #     namely the multiples of 1%·h^-1 (or 0.01 h^-1) within the range of 'zlim'
  #
  #   - the color scale 'Colors', comprising 'n.colors' colors. Setting 'colorPaletteOfGFtbox' to TRUE
  #     will use the color palette implemented in GFtbox. If FALSE, a blue-white-red colour scale will be used
  #     for negative-zero-postive values.

  if (missing(zlim))
    {
    if (round.zlim)
      {
      if (min(k, na.rm = T) < 0)
        {
        ifelse (Percent, mini <- floor(min(k, na.rm = T)*100)/100, mini <- floor(min(k, na.rm = T)))
        }
      else
        {
        mini <- round.value
        }
      if (max(k, na.rm = T) > 0)
        {
        ifelse (Percent, maxi <- ceiling(max(k, na.rm = T)*100)/100, maxi <- ceiling(max(k, na.rm = T)))
        }
      else
        {
        maxi <- round.value
        }
      }
    else
      {
      mini <- min(k, na.rm = T)
      maxi <- max(k, na.rm = T)
      }
    }
  else
    {
    mini <- zlim[1]
    maxi <- zlim[2]
    }
    
  if (Percent)
    {
    k <- k*100
    if (missing(zlim))
      {
      mini <- mini*100
      maxi <- maxi*100
      }
    }
    
  if (colorPaletteOfGFtbox)
    {  
    c(Colors, Range) := rainbowMap(c(mini, maxi), F, n.colors)
    }
  else
    {
    c(Colors, Range) := colorMap.Residual(c(mini, maxi), T, n.colors)
    }
  mini <- Range[1]
  maxi <- Range[2]

  below.mini <- which(k < mini)
  if (length(below.mini) > 0)
    {
    ifelse (fix.min, k[below.mini] <- mini, k[below.mini] <- NA)
    }
   
  above.maxi <- which(k > maxi)
  if (length(above.maxi) > 0)
    {
    ifelse (fix.max, k[above.maxi] <- maxi, k[above.maxi] <- NA)
    }

  if (mini != -Inf & mini != Inf & maxi != -Inf & maxi != Inf & maxi != mini)
    {
    ifelse (maxi - mini == 1, PowerOfTen <- 0, PowerOfTen <- ceiling(abs(log10(maxi - mini))) * log10(maxi - mini)/abs(log10(maxi - mini)))
    if (maxi - mini < 2 * 10^PowerOfTen) { PowerOfTen <- PowerOfTen - 1 }
    ticksGrowthScale <- seq(ceiling(mini/10^PowerOfTen)*10^PowerOfTen, floor(maxi/10^PowerOfTen)*10^PowerOfTen, 10^PowerOfTen)
    if (ticksGrowthScale[1] == mini) { ticksGrowthScale <- ticksGrowthScale[-1] }
    if (ticksGrowthScale[length(ticksGrowthScale)] == maxi) { ticksGrowthScale <- ticksGrowthScale[-length(ticksGrowthScale)] }
    }
  else
    {
    ticksGrowthScale <- NA
    }

  # --> edit 2015/12/16 #
  if (mini != 0 & mini == maxi)
    {
    Colors <- "#0000FF"
    }
  # edit 2015/12/16 <-- #

  return (list(k = k, zlim = c(mini, maxi), ticksGrowthScale = ticksGrowthScale, Colors = Colors))
  }


get.txt.legend <- function (k.growth, Percent)
  {
  # This function returns the label for the growth scale bar
  #
  # If 'Percent' is TRUE , a % character will be added in the label where relevant.
  # 'k.growth' is an integer between 2 and 11:
  #   - 2: Areal growth rate (k_area)
  #   - 3: Growth rate along the major axis (k_maj)
  #   - 4: Growth rate along the minor axis (k_min)
  #   - 5: Orientation of the major axis (theta)
  #   - 6: Rotation rate (phi)
  #   - 7: Anisotropy (1 - k_min/k_maj, dimensionless)
  #   - 8: Growth rate parallel to the leaf width which should have been made parallel to the X axis (kx, k_perpendicular_to_midline)
  #   - 9: Growth rate parallel to the midvein which should have been made parallel to the Y axis (ky, kml, k_midline)
  #   - 10: Cell area at time 1
  #   - 11: Cell area at time 2 (or sector area, as the cell may have divided)
 
  if (k.growth == 2)
    {
    ifelse (Percent, return (expression(paste("k"["area"], " (% h"^-1, ")"))),
                     return (expression(paste("k"["area"], " (h"^-1, ")"))))
    }
  else if (k.growth == 3)
    {
    ifelse (Percent, return (expression(paste("k"["maj"], " (% h"^-1, ")"))),
                     return (expression(paste("k"["maj"], " (h"^-1, ")"))))
    }
  else if (k.growth == 4)
    {
    ifelse (Percent, return (expression(paste("k"["min"], " (% h"^-1, ")"))),
                     return (expression(paste("k"["min"], " (h"^-1, ")"))))
    }
  else if (k.growth == 5)
    {
    return ("Orientation of the major axis (° to the x-axis)")
    }
  else if (k.growth == 6)
    {
    ifelse (Percent, return (expression(paste("Rotation rate (% h"^-1, ")"))),
                     return (expression(paste("Rotation rate (h"^-1, ")"))))
    }
  else if (k.growth == 7)
    {
    return ("Anisotropy")
    }
  else if (k.growth == 8)
    {
    ifelse (Percent, return (expression(paste("k"["per to midline"], " (% h"^-1, ")"))),
                     return (expression(paste("k"["per to midline"], " (h"^-1, ")"))))
    }
  else if (k.growth == 9)
    {
    ifelse (Percent, return (expression(paste("k"["midline"], " (% h"^-1, ")"))),
                     return (expression(paste("k"["midline"], " (h"^-1, ")"))))
    }
  else if (k.growth == 10)
    {
    return (expression(paste("Cell area (µm"^2, ")")))
    }
  else if (k.growth == 11)
    {
    return (expression(paste("Sector area (µm"^2, ")")))
    }
  }


get.txt.legend.Residual <- function (k.growth, Percent)
  {
  # Same as 'get.txt.legend()' for residual data

  if (k.growth == 2)
    {
    ifelse (Percent, return (expression(paste("Residual k"["area"], " (obs - fit, % h"^-1, ")"))),
                     return (expression(paste("Residual k"["area"], " (obs - fit, h"^-1, ")"))))
    }
  else if (k.growth == 3)
    {
    ifelse (Percent, return (expression(paste("Residual k"["maj"], " (obs - fit, % h"^-1, ")"))),
                     return (expression(paste("Residual k"["maj"], " (obs - fit, h"^-1, ")"))))
    }
  else if (k.growth == 4)
    {
    ifelse (Percent, return (expression(paste("Residual k"["min"], " (obs - fit, % h"^-1, ")"))),
                     return (expression(paste("Residual k"["min"], " (obs - fit, h"^-1, ")"))))
    }
  else if (k.growth == 5)
    {
    return ("Residual orientation of the major axis (obs - fit, ° to the x-axis)")
    }
  else if (k.growth == 6)
    {
    ifelse (Percent, return (expression(paste("Residual rotation rate (obs - fit, % h"^-1, ")"))),
                     return (expression(paste("Residual rotation rate (obs - fit, h"^-1, ")"))))
    }
  else if (k.growth == 7)
    {
    return ("Anisotropy")
    }
  else if (k.growth == 8)
    {
    ifelse (Percent, return (expression(paste("Residual k"["per to midline"], " (obs - fit, % h"^-1, ")"))),
                     return (expression(paste("Residual k"["per to midline"], " (obs - fit, h"^-1, ")"))))
    }
  else if (k.growth == 9)
    {
    ifelse (Percent, return (expression(paste("Residual k"["midline"], " (obs - fit, % h"^-1, ")"))),
                     return (expression(paste("Residual k"["midline"], " (obs - fit, h"^-1, ")"))))
    }
  else if (k.growth == 10)
    {
    return (expression(paste("Residual cell area (obs - fit, µm"^2, ")")))
    }
  else if (k.growth == 11)
    {
    return (expression(paste("Residual sector area (obs - fit, µm"^2, ")")))
    }
  }


plot.GrowthScale <- function (zlim, col.font, k.growth, Percent, txt.legend, Colors, drawTicks, ticksGrowthScale)
  {
  # This function plots the growth scale
  # and displays the legend 'txt.legend'
  
  y1 <- x1 <- 0
  y2 <- x2 <- 1
  xFactor <- 2
  par(mar = c(0,0,0,0))
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = c(x1 - x2*(xFactor-1), x2*xFactor), ylim = c(y1-0.5*(y2-y1), y2+1.75*(y2-y1)))
  if (missing(txt.legend)) { txt.legend <- get.txt.legend(k.growth, Percent) }
  slo <- (x2-x1)/length(Colors)
  int <- x1
  for (i in 0:(length(Colors)-1)) { rect(i*slo + int, y1, (i+1)*slo + int, y2, col = Colors[i+1], border = NA, xpd = NA) }
  rect(x1, y1, x2, y2, lwd = 2, xpd = NA, border = col.font)
  text((x1+x2)/2, y2+(y2-y1), txt.legend, xpd = NA, cex = 2, col = col.font)
  text(x1,(y1+y2)/2, round(zlim[1], 3), xpd = NA, cex = 2, pos = 2, col = col.font)
  text(x2,(y1+y2)/2, round(zlim[2], 3), xpd = NA, cex = 2, pos = 4, col = col.font)
  if (drawTicks & !is.na(ticksGrowthScale[1]))
    {
    ticksGrowthScale <- (ticksGrowthScale-zlim[1])/(zlim[2]-zlim[1])*length(Colors)*slo + int
    segments(ticksGrowthScale, y1, ticksGrowthScale, y1+0.1*(y2-y1), lwd = 2, col = col.font)
    segments(ticksGrowthScale, y2, ticksGrowthScale, y2-0.1*(y2-y1), lwd = 2, col = col.font)
    }
  }


# --> edit 2015/12/15 #
plot.DivisionScale <- function (zlim, col.font, txt.legend, Colors)
  {
  # This function plots the scale for the number of divisions (as integers)
  # and displays the legend 'txt.legend'

  y1 <- x1 <- 0
  y2 <- x2 <- 1
  xFactor <- 1.5
  par(mar = c(0,0,0,0))
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = c(x1 - x2*(xFactor-1), x2*xFactor), ylim = c(y1-0.5*(y2-y1), y2+1.75*(y2-y1)))
  boxes <- (zlim[1]):(zlim[2])
  slo <- (x2-x1)/length(boxes)
  int <- x1
  for (i in 0:(length(boxes)-1))
    {
    ifelse (zlim[1] != 0 & zlim[1] == zlim[2],
            color <- Colors,
            color <- Colors[round((1 - zlim[1]*(length(Colors)-1)/(zlim[2]-zlim[1])) + (length(Colors)-1)/(zlim[2]-zlim[1]) * boxes[i+1])])
    rect(i*slo + int, y1, (i+1)*slo + int, y2, col = color, lwd = 2, xpd = NA, border = col.font)
    text((i+0.5)*slo + int,(y1+y2)/2, boxes[i+1], xpd = NA, cex = 1.7, col = col.font)
    }
  text((x1+x2)/2, y2+(y2-y1), txt.legend, xpd = NA, cex = 2, col = col.font)
  }
# edit 2015/12/15 <-- #


plot.DivAndCompLegend <- function (div.parameter, col.font, ClassName)
  {
  # This function plots the legend
  # for cell division and/or competence.
  # Only 1 class is currently allowed.
  
  y1 <- x1 <- 0
  y2 <- x2 <- 1
  xFactor <- 2
  xmin <- x1 - x2*(xFactor-1)
  xmax <- x2*xFactor
  ymin <- y1-0.5*(y2-y1)
  ymax <- y2+1.75*(y2-y1)
  par(mar = c(0,0,0,0), xpd = NA)
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = c(xmin, xmax), ylim = c(ymin, ymax))

  y1 <- 1.9
  y2 <- 1.2
  txt <- c("Cell division", "Competence for division")
  col <- c("olivedrab2", "olivedrab") # colors for competence and division have been switched on 2016/01/07 #
  if (div.parameter %in% c("div", "comp"))
    {
    if (is.null(ClassName))
      {
      ifelse (div.parameter == "div", idx <- 1, idx <- 2)
      ifelse (div.parameter == "div", x1 <- 0.05, x1 <- -0.35)
      rect(x1, y1, x1+0.1, y2, border = "grey", col = col[idx])
      text(x1+0.1, (y1+y2)/2, txt[idx], cex = 2, pos = 4, col = col.font)
      }
    else
      {
      xseq <- seq(xmin + 0.25*abs(xmin), xmax + 0.25*abs(xmin), length.out = nrow(ClassName)+1)
      xseq <- xseq[1:nrow(ClassName)] + diff(xseq)/2
      yseq <- seq(ymin + 0.25*abs(ymin), ymax - 0.2*abs(ymax), length.out = 3)
      yseq <- yseq[1:2] + diff(yseq)/2
      text(xseq, ymax - 0.1*abs(ymax), ClassName$ClassName, cex = 1, col = col.font)
      text(xmin - 0.15 * abs(xmin), yseq[2], ifelse (div.parameter == "div", "Has divided", "Is competent"), cex = 1, pos = 4, col = col.font)
      ifelse (div.parameter == "div", color <- ClassName$DivisionColor, color <- ClassName$CompetenceColor)
      rect(xseq-0.03*(xmax-xmin), yseq[2]-0.15*yseq[2], xseq+0.03*(xmax-xmin), yseq[2]+0.15*yseq[2], border = "grey", col = color)
      }
    }
  else
    {
    if (is.null(ClassName))
      {
      x1 <- c(-1, 0.15)
      rect(x1, y1, x1+0.1, y2, border = "grey", col = col)
      text(x1+0.1, (y1+y2)/2, txt, cex = 2, pos = 4, col = col.font)
      }
    else
      {
      xseq <- seq(xmin + 0.25*abs(xmin), xmax + 0.25*abs(xmin), length.out = nrow(ClassName)+1)
      xseq <- xseq[1:nrow(ClassName)] + diff(xseq)/2
      yseq <- seq(ymin + 0.25*abs(ymin), ymax - 0.2*abs(ymax), length.out = 3)
      yseq <- yseq[1:2] + diff(yseq)/2
      text(xseq, ymax - 0.1*abs(ymax), ClassName$ClassName, cex = 1, col = col.font)
      text(xmin - 0.15 * abs(xmin), yseq, c("Is competent", "Has divided"), cex = 1, pos = 4, col = col.font)
      rect(xseq-0.03*(xmax-xmin), yseq[2]-0.15*yseq[2], xseq+0.03*(xmax-xmin), yseq[2]+0.15*yseq[2], border = "grey", col = ClassName$DivisionColor)
      rect(xseq-0.03*(xmax-xmin), yseq[1]-0.15*yseq[2], xseq+0.03*(xmax-xmin), yseq[1]+0.15*yseq[2], border = "grey", col = ClassName$CompetenceColor)
      }
    }    
  }


plot.ScaleBar <- function (ylim, scaleBar, tick, col.font)
  {
  # This function plots the scale bar in the lowest margin
  # given as a numeric value in as 'scaleBar'.
  # The argument 'ylim' should be in the form:
  #     ylim = c(lowest_value, highest_value)
  
  yScaleBar <- ylim[1] - (ylim[2] - ylim[1]) * 0.07     # position of the bar
  yScaleValue <- ylim[1] - (ylim[2] - ylim[1]) * 0.125  # position of the value
  arrows(+0.1*scaleBar, yScaleBar, -scaleBar/2, yScaleBar, angle = 90, lwd = 2, length = ifelse (tick, 0.05, 0), xpd = NA, col = col.font)
  arrows(-0.1*scaleBar, yScaleBar, +scaleBar/2, yScaleBar, angle = 90, lwd = 2, length = ifelse (tick, 0.05, 0), xpd = NA, col = col.font)
  text(0, yScaleValue, paste(scaleBar, "µm"), xpd = NA, cex = 1.5, col = col.font)
  }


plot.Anisotropy <- function (aniso.threshold, aniso.lwd, aniso.lwd.constant, AreaTracked, Shapes, ShapesX, ShapesY, Growth, meanAngle)
  {
  # This function plots anisotropy over the surface of the already plotted cells.
  # Anisotropy will be plotted only where superior to 'aniso.threshold',
  # which should lie between 0 and 1.
  #
  # Regarding the length and thickness of the lines:
  #
  #   - The length of each line is automatically calculated as a function
  #     of the length of each cell major axis of growth.
  #
  #   - Except if 'aniso.lwd.constant' is TRUE, the thickness of each line is directly proportional
  #     to the anisotropy value. Unless provided in 'aniso.lwd', the factor by which the anisotropy
  #     value will be multiplied is tentatively calculated using an empirical function of
  #       i) the number of cells tracked 'nlevels(Shapes$Cell)' and
  #       ii) the whole area tracked 'AreaTracked'.

  if (missing(aniso.lwd)) { aniso.lwd <- max(3, 9 - 8e-07 * AreaTracked * nlevels(Shapes$Cell)) }
  for (cell in levels(Shapes$Cell))
    {
    if (Growth[Growth$Cell == cell, 7] > aniso.threshold)
      {
      Xs <- as.numeric(ShapesX[ShapesX$Cell == cell, 2:ncol(ShapesX)])
      Ys <- as.numeric(ShapesY[ShapesY$Cell == cell, 2:ncol(ShapesY)])
      c(Xm, Ym) := get.Centroid(cbind(na.omit(Xs), na.omit(Ys)))
      Angle <- -Growth[Growth$Cell == cell, 5] # Warning: the angle in the output growth file is negatively oriented #
      Angle <- Angle + meanAngle
      Ys.rot <- -(Xm - Xs) * sin(Angle*pi/180) + (Ym - Ys) * cos(Angle*pi/180)
      r <- 0.6 * sqrt((Xs[which.min(Ys.rot^2)] - Xm)^2 + (Ys[which.min(Ys.rot^2)] - Ym)^2)
      x1 <- Xm + r*cos(Angle*pi/180)
      x2 <- Xm - r*cos(Angle*pi/180)
      y1 <- Ym + r*sin(Angle*pi/180)
      y2 <- Ym - r*sin(Angle*pi/180)
      ifelse (aniso.lwd.constant, lwd <- aniso.lwd * 0.5, lwd <- aniso.lwd * Growth[Growth$Cell == cell, 7])
      segments(x1, y1, x2, y2, lwd = lwd)
      }
    }
  }


plot.growth <- function(InfoVertices, Vertices, Cells, Divisions, meanAngle, Growth, Shapes, 
                        alignToPetioleLaminaBoundary = T, cellAtPetioleLaminaBoundary,
                        Image1, Image2, leafShape, before = F, k.growth = 2,
                        PNG = F, black = T, xlim, ylim, scaleBar, tick = F,
                        Percent, round.zlim = T, zlim, fix.min = T, fix.max = T,
                        colorPaletteOfGFtbox = T, growthScale = T, drawTicks = T, txt.legend,
                        anisotropy, aniso.threshold = 0.05, aniso.lwd, aniso.lwd.constant = F)
  {
  # This function plots growth between 'Image 1' and 'Image2'
  # either computed in PointTracker or recomputed in R.
  #
  # The following arguments are from 'import.processedData()' (if processed in PointTracker),
  # or 'from import.rawData()' followed by 'recompute.GrowthAndShapes()' (if recomputed in R):
  #
  #     'InfoVertices', 'Vertices', 'Cells', 'Divisions', 'meanAngle', 'Growth', 'Shapes'
  #
  # Other parametres may be modified by the user:
  #
  #   - Set 'alignToPetioleLaminaBoundary' to TRUE to offset the coordinates to the petiole-lamina boundary,
  #     possibly entered as a cell number using 'cellAtPetioleLaminaBoundary'.
  #     See 'alignToPLB()' and 'scale.Objects()' for more details.
  #
  #   - Setting 'leafShape' to TRUE will draw the leaf outline
  #     saved as XY coordinates using ImageJ under 'Leaf Shape\Name_Of_Image_i.txt'.
  #
  #   - Setting 'before' to FALSE (default) will draw cells
  #    (and leaf outline if 'leafShape' is TRUE) using the second image ('Image2').
  #    If set to FALSE, the first image ('Image1') is used.
  #
  #   - 'k.growth' is an integer between 2 and 11 which gives the growth variable
  #     to be plotted as explained in 'get.txt.legend()'. The default (2) is 'k_area'.
  #
  #   - Setting 'PNG' to TRUE will prevent to open an on-screen device,
  #     and thus allow to plot the figure in an alternative graphics device such as 'png()'.
  #
  #   - Setting 'black' to TRUE (default) will draw the leaf on a black background.
  #     White background if FALSE.
  #
  #   - 'xlim' and 'ylim' are the limits for the x and y axes of the main plot, respectively.
  #     There should be either missing (and their computation will be automatic)
  #     or given in microns relative to the petiole lamina-boundary as vectors
  #     with two numeric values, namely:
  #         'xlim = c(x.min, x.max)'
  #         'ylim = c(y.min, y.max)'
  #     To obtain plots of several images at the same scale, run the function for each image
  #     by setting a constant 'ylim' to the optimal boundaries of the longest leaf, e.g. 'ylim = c(-200, 1000)'.
  #     In all cases the X/Y ratio will be constrained to 1, so defining 'xlim' is generally unuseful.
  #
  #   - 'scaleBar' is a numeric value in microns giving the length of the scale bar.
  #     For example, use 'scaleBar = 100' to draw a scale bar of 100 microns.
  #     If missing, the scale bar is not drawn.
  #
  #   - Setting 'tick' to FALSE (default) will make the scale bar as a a simple segment.
  #     Setting it to TRUE will draw the ticks on each side of the scale bar.
  #
  #   - Setting 'Percent' to TRUE will make growth to be expressed in %·h^-1.
  #     If missing, automatically set to TRUE (except if the variable plotted
  #     is anisotropy, the orientation of the major axis, or cell area).
  #
  #   - 'round.zlim', 'zlim', 'fix.min', 'fix.max': see 'set.zlim()'.
  #
  #   - Setting 'colorPaletteOfGFtbox' to TRUE (the default) will use the color palette
  #     implemented in GFtbox.
  #
  #   - Setting 'growthScale' to TRUE (the default) will plot the growth scale using 'plot.GrowthScale()'.
  #
  #   - Setting 'drawTicks' to TRUE (the default) will draw ticks on the growth scale
  #     at the values set in 'ticksGrowthScale' (see 'set.zlim()').
  #
  #   - 'txt.legend' is the text above the growth scale. Can be of type 'character'
  #     or 'expression' (the latter being usefull for subscripts, superscripts, greek letters...).
  #     If missing, it is automatically obtained using 'get.txt.legend()'.
  #
  #   - Setting 'anisotropy' to TRUE will make anisotropy to be plotted using 'plot.Anisotropy()'.
  #     If missing, automatically set to TRUE only if 'k_area' is plotted.
  #
  #   - 'aniso.threshold' (default = 0.05 = 5%), 'aniso.lwd', 'aniso.lwd.constant' (default = F):
  #     see 'plot.Anisotropy()'.


  # Perform initial checks on parameters given by the user
  initialChecks(scaleBar, xlim, ylim, zlim, aniso.threshold, aniso.lwd, check.type = "growth")

  # Prepare the plot window
  if (!PNG) { windows() }
  if (black) { par(bg = "black") }
  ifelse (black, col.font <- "white", col.font <- "black")
  layout(matrix(c(2, 1)), heights = c(1, 5))
  par(mar = c(4,0,0,0))
  
  # Import leaf outline, cell shapes and vertices
  if (missing(leafShape))
    {
    ifelse (file.exists(paste("Leaf Shape/", substr(Image1, 1, nchar(Image1) - 3), "txt", sep = "")) &
            file.exists(paste("Leaf Shape/", substr(Image2, 1, nchar(Image2) - 3), "txt", sep = "")),
            leafShape <- T, leafShape <- F)
    }
  c(VertX, VertY, ShapesX, ShapesY, LeafShape, LeafSpline) := scale.Objects(Image1, Image2, before, meanAngle, leafShape, alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary, Shapes, Cells, Vertices, InfoVertices)
  
  # Leaf outline
  if (leafShape)
    {
    if (missing(xlim)) { xlim <- c(-1.01*max(abs(LeafSpline$x)), 1.01*max(abs(LeafSpline$x))) }
    if (missing(ylim)) { ylim <- c(min(LeafSpline$y) - 0.01*abs(min(LeafSpline$y)), 1.01*max(LeafSpline$y)) }
    image(x = 0, y = 0, z = matrix(1,1), col = NA, asp = 1, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = xlim, ylim = ylim) # image() is used so that axes limits match those for kriging
    #plot(VertX, VertY, asp = 1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = ylim)
    polygon(LeafSpline, lwd = 2, border = "grey")
    }
  else
    {
    if (missing(xlim)) { xlim <- c(min(VertX, na.rm = T), max(VertX, na.rm = T)) }
    if (missing(ylim)) { ylim <- c(min(VertY, na.rm = T) - 0.01*abs(min(VertY, na.rm = T)), 1.01*max(VertY, na.rm = T)) }
    image(x = 0, y = 0, z = matrix(1,1), col = NA, asp = 1, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = xlim, ylim = ylim) # image() is used so that axes limits match those for kriging
    #plot(VertX, VertY, asp = 1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    }

  # Growth
  k <- Growth[, k.growth]
  if (k.growth == 5) { k <- -k + meanAngle }
  if (missing(Percent)) { ifelse (k.growth %in% c(5, 7, 10, 11), Percent <- F, Percent <- T) }
  c(k, zlim, ticksGrowthScale, Colors) := set.zlim(k, Percent, round.zlim, zlim, fix.min, fix.max, colorPaletteOfGFtbox)
  mini <- zlim[1]
  maxi <- zlim[2]
      
  # Cells
  AreaTracked <- 0
  for (cell in levels(Shapes$Cell))
    {
    color <- Colors[round((1 - mini*(length(Colors)-1)/(maxi-mini)) + (length(Colors)-1)/(maxi-mini) * k[Growth$Cell == cell])]
    vert <- na.omit(as.numeric(Cells[Cells$Cell == cell, -1]))
    XY <- cell.coord(vert, VertX, VertY, Vertices)
    polygon(XY, col = color, border = "grey")
    AreaTracked <- AreaTracked + areapl(na.omit(XY))
    }

  # Anisotropy
  if (missing(anisotropy)) { ifelse (k.growth == 2, anisotropy <- T, anisotropy <- F) }
  if (anisotropy)
    {
    plot.Anisotropy(aniso.threshold, aniso.lwd, aniso.lwd.constant, AreaTracked, Shapes, ShapesX, ShapesY, Growth, meanAngle)
    }
    
  # Scale bar
  if (!missing(scaleBar)) { plot.ScaleBar(ylim, scaleBar, tick, col.font) }
  
  # Growth scale
  if (growthScale) { plot.GrowthScale(zlim, col.font, k.growth, Percent, txt.legend, Colors, drawTicks, ticksGrowthScale) }
  }


plot.kriging <- function(InfoVertices, Vertices, Cells, Divisions, meanAngle, Growth, Shapes, 
                         alignToPetioleLaminaBoundary = T, cellAtPetioleLaminaBoundary,
                         Image1, Image2, leafShape, leafShapeKriging = F, before = F, k.growth = 2,
                         PNG = F, black = T, xlim, ylim, scaleBar, tick = F,
                         Percent, round.zlim = T, zlim, fix.min = T, fix.max = T,
                         colorPaletteOfGFtbox = T, growthScale = T, drawTicks = T, txt.legend,
                         anisotropy = F, aniso.threshold = 0.05, aniso.lwd, aniso.lwd.constant = F,
                         n.pred = 100, polynomDegree = 3, rangeParameter = 0.001, contour.levels,
                         cellularScale = F, plotResidual = F, exportCellValues = F)
  {
  # This function performs and plots a kriging estimation of growth
  # between 'Image 1' and 'Image2' either computed in PointTracker or recomputed in R.
  #
  # Arguments are the same as for 'plot.growth()', plus some extra:
  #
  #   - Setting 'leafShapeKriging' to TRUE will perform a kriging prediction over the whole leaf surface,
  #     whichs yields generally a bad result except if the major part of the leaf was tracked.
  #     If FALSE (default), the prediction will be restricted to the tracked area.
  #
  #   - The default 'anisotropy' is always set to FALSE, because it is difficult to see both the anisotropy lines
  #     and the contours. However, it still can be plotted as in 'plot.growth()' by setting 'anisotropy' to TRUE.
  #
  #   - 'n.pred' is the number of x-coordinates and y-coordinates where the kriging prediction is performed
  #     over the specified surface. The total number of predictions is n.pred^2.
  #     'n.pred' should be an integer between 10 (very pixelated) and 1000 (very smooth).
  #     Default is 100, as a good compromise between quality of the prediction and time required
  #     to perform and plot the estimated surface.
  #
  #   - 'polynomDegree' is the degree of polynomial surface (see '?surf.gls'). It should be an integer between 2 and 6.
  #     Increasing 'polynomDegree' will yield more convoluted profiles. Setting it to maximum does not generally
  #     produce graphs making sense. Default is 3.
  #
  #   - 'rangeParameter' is the range used in the spatial exponental covariance function (see '?expcov').
  #     It should be a strictly posititve number. Default is 0.001.
  #     Increasing 'rangeParameter' will yield more convoluted profiles.
  #
  #   - 'contour.levels' is a numeric vector of at least one value giving the levels where the contours are plotted.
  #     There should be given in %·h^-1 if Percent is TRUE, in h^-1 otherwise. Examples are: 
  #       - 'contour.levels = 0.15'
  #       - 'contour.levels = c(0, 1, 4)'
  #       - 'contour.levels = seq(0, 0.06, 0.01)'
  #     If missing, all multiples of 1 %·h^-1 (or 0.01 h^-1) within the range of 'zlim' will be plotted.
  #
  #   - Setting 'cellularScale' to TRUE will produce a map
  #     of the kriging estimate averaged over the area of each cell.
  #
  #   - Setting 'plotResidual' to TRUE will produce a map
  #     of the difference between the measured and estimated value for each cell (residuals).
  #
  #   - Setting 'exportCellValues' to TRUE will return a dataframe
  #     of the cell measured and actual values, as well as their difference (residuals).

  # Perform initial checks on parameters given by the user
  initialChecks(scaleBar, xlim, ylim, zlim, aniso.threshold, aniso.lwd,
                n.pred, polynomDegree, rangeParameter, contour.levels,
                leafShapeKriging, cellularScale, plotResidual, exportCellValues,
                check.type = "kriging")

  # Prepare the plot window
  if (!PNG) { windows() }
  if (black) { par(bg = "black") }
  ifelse (black, col.font <- "white", col.font <- "black")
  layout(matrix(c(2, 1)), heights = c(1, 5))
  par(mar = c(4,0,0,0))
  
  # Import leaf outline, cell shapes and vertices
  if (missing(leafShape))
    {
    ifelse (file.exists(paste("Leaf Shape/", substr(Image1, 1, nchar(Image1) - 3), "txt", sep = "")) &
            file.exists(paste("Leaf Shape/", substr(Image2, 1, nchar(Image2) - 3), "txt", sep = "")),
            leafShape <- T, leafShape <- F)
    }
  c(VertX, VertY, ShapesX, ShapesY, LeafShape, LeafSpline) := scale.Objects(Image1, Image2, before, meanAngle, leafShape, alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary, Shapes, Cells, Vertices, InfoVertices)

  # Associate centroid (Xm, Ym) to growth
  centro <- data.frame(x = NULL, y = NULL, z = NULL)
  for (cell in levels(Shapes$Cell))
    {
    vert <- na.omit(as.numeric(Cells[Cells$Cell == cell, -1]))
    XY <- cell.coord(vert, VertX, VertY, Vertices)
    c(Xm, Ym) := get.Centroid(XY)
    K <- Growth[Growth$Cell == cell, k.growth]
    if (k.growth == 5) { K <- -K + meanAngle }
    centro <- rbind(centro, data.frame(x = Xm, y = Ym, z = K))
    }

  # Kriging over the whole leaf or the tracked region
  ifelse (leafShapeKriging, minX <- min(LeafSpline$x), minX <- min(VertX, na.rm = T))
  ifelse (leafShapeKriging, maxX <- max(LeafSpline$x), maxX <- max(VertX, na.rm = T))
  ifelse (leafShapeKriging, minY <- min(LeafSpline$y), minY <- min(VertY, na.rm = T))
  ifelse (leafShapeKriging, maxY <- max(LeafSpline$y), maxY <- max(VertY, na.rm = T))
  kr <- surf.gls(polynomDegree, expcov, centro, d = rangeParameter)
  prsurf <- prmat(kr, minX, maxX, minY, maxY, n.pred-1)

  # Find which predicted values fall in the whole leaf or the tracked region
  matpred <- expand.grid(y = prsurf$y, x = prsurf$x, KEEP.OUT.ATTRS = F)
  if (leafShapeKriging)
    {
    pip <- matrix(point.in.polygon(matpred$x, matpred$y, LeafSpline$x, LeafSpline$y), nrow = n.pred, byrow = T)
    }
  else
    {
    pip <- matrix(nrow = n.pred, ncol = n.pred)
    if (cellularScale | plotResidual) { kCellSmooth <- vector() }
    for (cell in levels(Shapes$Cell))
      {
      vert <- na.omit(as.numeric(Cells[Cells$Cell == cell, -1]))
      XY <- cell.coord(vert, VertX, VertY, Vertices)
      pip.cell <- pip.temp <- matrix(point.in.polygon(matpred$x, matpred$y, XY[,1], XY[,2]), nrow = n.pred, byrow = T)
      if (cellularScale | plotResidual) { kCellSmooth <- c(kCellSmooth, mean(prsurf$z[pip.cell > 0])) }
      pip.temp[pip.temp > 0] <- 1
      pip.temp[pip == 1] <- 1
      pip <- pip.temp
      }
    }
  prsurf$z[pip == 0] <- NA

  # Plot the predicted surface
  if (missing(Percent)) { ifelse (k.growth %in% c(5, 7, 10, 11), Percent <- F, Percent <- T) }

  if (cellularScale | plotResidual)
    {
    if (!plotResidual)
      {
      c(k, zlim, ticksGrowthScale, Colors) := set.zlim(kCellSmooth, Percent, round.zlim, zlim, fix.min, fix.max, colorPaletteOfGFtbox)
      c(prsurf$z, zlim, ticksGrowthScale, Colors) := set.zlim(k = prsurf$z, Percent, round.zlim, zlim, fix.min, fix.max, colorPaletteOfGFtbox) # for contours
      }
    else
      {
      c(k, zlim, ticksGrowthScale, Colors) := set.zlim(centro$z - kCellSmooth, Percent, round.zlim = T, zlim, fix.min, fix.max, colorPaletteOfGFtbox = F)
      }
    mini <- zlim[1]
    maxi <- zlim[2]
    COLORS <- NA
    }
  else
    {
    c(prsurf$z, zlim, ticksGrowthScale, Colors) := set.zlim(k = prsurf$z, Percent, round.zlim, zlim, fix.min, fix.max, colorPaletteOfGFtbox)
    COLORS <- Colors
    color <- NA
    }

  if (leafShape | leafShapeKriging)
    {
    if (missing(xlim)) { xlim <- c(-1.01*max(abs(LeafSpline$x)), 1.01*max(abs(LeafSpline$x))) }
    if (missing(ylim)) { ylim <- c(min(LeafSpline$y) - 0.01*abs(min(LeafSpline$y)), 1.01*max(LeafSpline$y)) }
    image(prsurf, col = COLORS, asp = 1, bty = "n", xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, zlim = zlim)
    if (leafShape) { polygon(LeafSpline, lwd = 2, border = "grey") }
    }
  else
    {
    if (missing(xlim)) { xlim <- c(min(VertX, na.rm = T), max(VertX, na.rm = T)) }
    if (missing(ylim)) { ylim <- c(min(VertY, na.rm = T) - 0.01*abs(min(VertY, na.rm = T)), 1.01*max(VertY, na.rm = T)) }
    image(prsurf, col = COLORS, asp = 1, bty = "n", xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, zlim = zlim)
    }
  
  # Cell outlines
  AreaTracked <- 0
  for (cell in levels(Shapes$Cell))
    {
    if (cellularScale | plotResidual) { color <- Colors[round((1 - mini*(length(Colors)-1)/(maxi-mini)) + (length(Colors)-1)/(maxi-mini) * k[Growth$Cell == cell])] }
    vert <- na.omit(as.numeric(Cells[Cells$Cell == cell, -1]))
    XY <- cell.coord(vert, VertX, VertY, Vertices)
    polygon(XY, border = "grey", col = color)
    AreaTracked <- AreaTracked + areapl(na.omit(XY))
    }
  
  # Contours
  if (!plotResidual)
    {
    if (missing(contour.levels))
      {
      ifelse (Percent, contour.levels <- seq(ceiling(zlim[1]), floor(zlim[2]), 1), contour.levels <- seq(ceiling(zlim[1]*100)/100, floor(zlim[2]*100)/100, 0.01))
      }
    contour.levels <- unique(contour.levels)
    contour(prsurf, levels = contour.levels, add = T, lwd = 2)
    }
  
  # Anisotropy
  if (anisotropy)
    {
    plot.Anisotropy(aniso.threshold, aniso.lwd, aniso.lwd.constant, AreaTracked, Shapes, ShapesX, ShapesY, Growth, meanAngle)
    }
    
  # Scale bar
  if (!missing(scaleBar)) { plot.ScaleBar(ylim, scaleBar, tick, col.font) }
  
  # Growth scale
  if (missing(txt.legend)) { txt.legend <- get.txt.legend(k.growth, Percent) }
  if (plotResidual) { txt.legend <-  get.txt.legend.Residual(k.growth, Percent) }
  if (growthScale) { plot.GrowthScale(zlim, col.font, k.growth, Percent, txt.legend, Colors, drawTicks, ticksGrowthScale) }

  # Data export
  if ((cellularScale | plotResidual) & exportCellValues) { return (CellValues = data.frame(Cell = levels(Shapes$Cell), K = K, kCellSmooth = kCellSmooth, Residual = K - kCellSmooth)) }

  }


plot.division <- function(InfoVertices, Vertices, Cells, Divisions, meanAngle, Shapes, Div,
                          alignToPetioleLaminaBoundary = T, cellAtPetioleLaminaBoundary,
                          Image1, Image2, leafShape, before = F,
                          PNG = F, black = T, xlim, ylim, scaleBar, tick = F,
                          div.parameter = "div&comp", round.zlim = T, zlim, fix.min = T, fix.max = T,
                          #colorPaletteOfGFtbox = T, growthScale = T, drawTicks = T, txt.legend)
                          colorPaletteOfGFtbox = T, growthScale = T, drawTicks = T, txt.legend, ini = 1, show.cell.number = F) # edits 2014/09/11 and 2016/01/15
  {
  # This function plots the division and competence class of cells computed
  # between 'Image 1' and 'Image2' using 'process.DivWithInterval()' which produces
  # the argument 'Div'. All other arguments are the same as for 'plot.growth()'.
  #
  # Furthermore, 'div.parameter' provides the division parameter
  # to be plotted and should be given as one of the following character strings:
  #
  #   - "div&comp" (the default): both division and competence (if a cell has divided
  #     and is competent for division, the color for division will dominate)
  #
  #   - "div": division only (if a cell has divided, it will have a color)
  #
  #   - "comp": competence only (if a cell is competent for division, it will have a color)
  #
  #   - "progeny": each cell lineage/clone/sector/progeny will be shown
  #
  #   - "CellArea": area of the cell (at 'Image1' if 'before' is TRUE, at 'Image2' if 'before' is FALSE).
  #     The area of every cells will be plotted.
  #
  #   - "AreaAtTime1": area of the cell at 'Image1'.
  #     It makes more sense to use it rather than "AreaAtTime2" if 'before' is TRUE.
  #     The area of every cells will be plotted.
  #
  #   - "AreaAtTime2": area of the cell (or daughter sector) at 'Image2'.
  #     It makes more sense to use it rather than "AreaAtTime1" if 'before' is FALSE.
  #     The area of every sectors will be plotted.
  #
  #   - "Area_At_Division": the area of the cell (if 'before' is TRUE) or of the mother cell (if 'before' is FALSE)
  #     at the time of division is displayed according to a color scale
  #
  #   - "Symmetry": the symmetry of division (computed as the minimal area of the two daughter cells
  #     divided by half of the area of the mother, between 0 and 1) is displayed according to a color scale
  #
  #   - "Cell_Cycle_Duration": the duration of cell cycle (if date of birth available)
  #     is displayed according to a color scale
  #
  #   - "Number_Of_Divisions": the number of divisions for each cell between 'Image1' and 'Image2'.
  #     If 'before' is FALSE, daughter cells are plotted with the number of cell division that actually relates to each cell.
  #     If 'before' is TRUE, mother cells are plotted with the maximal number of cell division is found for each lineage.
  #
  #   - 'show.cell.number' will display the number of each cell (as in PointTracker) in the font color if set to TRUE
  #
  # Note that the color scale will be plotted only if 'div.parameter' is either
  # "AreaAtTime1", "AreaAtTime2", "Area_At_Division", "Symmetry" or "Cell_Cycle_Duration".
  # This means that the arguments 'round.zlim', 'zlim', 'fix.min', 'fix.max', and 'txt.legend' will be used
  # only in these five cases. Otherwise, a legend will be displayed for division and/or competence (nothing for progeny).
  #
  # Note also that 'ini' is used only for progeny display (see 'process.AllDiv()').

  # Perform initial checks on parameters given by the user
  initialChecks(scaleBar, xlim, ylim, zlim, div.parameter = div.parameter, check.type = "division")

  # Prepare the plot window
  if (!PNG) { windows() }
  if (black) { par(bg = "black") }
  ifelse (black, col.font <- "white", col.font <- "black")
  layout(matrix(c(2, 1)), heights = c(1, 5))
  par(mar = c(4,0,0,0))
  
  # Import leaf outline, cell shapes and vertices
  if (missing(leafShape))
    {
    ifelse (file.exists(paste("Leaf Shape/", substr(Image1, 1, nchar(Image1) - 3), "txt", sep = "")) &
            file.exists(paste("Leaf Shape/", substr(Image2, 1, nchar(Image2) - 3), "txt", sep = "")),
            leafShape <- T, leafShape <- F)
    }
  c(VertX, VertY, ShapesX, ShapesY, LeafShape, LeafSpline) := scale.Objects(Image1, Image2, before, meanAngle, leafShape, alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary, Shapes, Cells, Vertices, InfoVertices)
  
  # Leaf outline
  if (leafShape)
    {
    if (missing(xlim)) { xlim <- c(-1.01*max(abs(LeafSpline$x)), 1.01*max(abs(LeafSpline$x))) }
    if (missing(ylim)) { ylim <- c(min(LeafSpline$y) - 0.01*abs(min(LeafSpline$y)), 1.01*max(LeafSpline$y)) }
    image(x = 0, y = 0, z = matrix(1,1), col = NA, asp = 1, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = xlim, ylim = ylim) # image() is used so that axes limits match those for kriging
    #plot(VertX, VertY, asp = 1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = ylim)
    polygon(LeafSpline, lwd = 2, border = "grey")
    }
  else
    {
    if (missing(xlim)) { xlim <- c(min(VertX, na.rm = T), max(VertX, na.rm = T)) }
    if (missing(ylim)) { ylim <- c(min(VertY, na.rm = T) - 0.01*abs(min(VertY, na.rm = T)), 1.01*max(VertY, na.rm = T)) }
    image(x = 0, y = 0, z = matrix(1,1), col = NA, asp = 1, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = xlim, ylim = ylim) # image() is used so that axes limits match those for kriging
    #plot(VertX, VertY, asp = 1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    }

  # Set limits of the division parameter to be plotted (if required)
  if (div.parameter %in% c("CellArea", "AreaAtTime1", "AreaAtTime2", "Area_At_Division", "Symmetry", "Cell_Cycle_Duration", "Number_Of_Divisions")) # edit 2015/12/15
    {
    if (div.parameter == "CellArea")
      {
      d <- vector()
      for (i in 1:nrow(Div))
        {
        cell <- as.character(Div$Cell[i])
        vert <- na.omit(as.numeric(Cells[Cells$Cell == cell, -1]))
        XY <- cell.coord(vert, VertX, VertY, Vertices)
        d <- c(d, areapl(na.omit(XY)))
        }
      }
    else
      {
      d <- Div[, div.parameter]
      }
    # --> edit 2015/12/15 #
    if (div.parameter == "Number_Of_Divisions")
      {
      if (missing(zlim))
        {
        round.zlim <- T
        ifelse (length(which(!is.na(d))) == 0, round.value <- 0, # if no cell has divided, only zero will be plotted on the legend
                                               round.value <- 1) # if at least one cell has divided, the default minimum will be one
        }
      else
        {
        zlim[1] <- max(c(1, round(zlim[1])))
        zlim[2] <- round(zlim[2])
        #fix.min <- F # otherwise, if zlim[1] > 1, it doesn't make sense for cells that divide less than zlim[1] but more than zero
        }
      }
    else
      {
      round.value <- 0
      }
    c(d, zlim, ticksGrowthScale, Colors) := set.zlim(d, Percent = F, round.zlim, zlim, fix.min, fix.max, colorPaletteOfGFtbox, round.value = round.value)
    #c(d, zlim, ticksGrowthScale, Colors) := set.zlim(d, Percent = F, round.zlim, zlim, fix.min, fix.max, colorPaletteOfGFtbox)
    # edit 2015/12/15 <-- #
    mini <- zlim[1]
    maxi <- zlim[2]
    }

  # Load user's colors (if any)
  else if (div.parameter %in% c("div&comp", "div", "comp"))
    {
    ClassName <- NULL
    if (file.exists("Divisions/DivisionClass.xlsx"))
      {
      if (!"XLConnect" %in% installed.packages()[, "Package"]) { install.packages("XLConnect") }
      if (!"package:XLConnect" %in% search()) { library("XLConnect", character.only = T) } # the package 'XLConnect' cannot be loaded at the beginning as it impairs 'rJava'
      workBook <- loadWorkbook("Divisions/DivisionClass.xlsx")
      sheetNames <- getSheets(workBook)
      if ("DivClass" %in% sheetNames)
        {
        DivType <- readWorksheet(workBook, sheet = "DivClass")
        if ("ClassName" %in% sheetNames)
          {
          ClassName <- readWorksheet(workBook, sheet = "ClassName")
          if (nlevels(as.factor(DivType$Class)) > nlevels(as.factor(ClassName$Class)))
            {
            stop("The classes in the sheet 'ClassName' do not correspond to the classes in the sheet 'DivClass'")
            }
          else
            {
            for (i in 1:nlevels(as.factor(DivType$Class)))
              {
              matchClass <- which(levels(as.factor(ClassName$Class)) == levels(as.factor(DivType$Class))[i])
              if (length(matchClass) != 1) { stop("The classes in the sheet 'ClassName' do not correspond to the classes in the sheet 'DivClass'") }
              }
            }
          }
        else
          {
          ClassName <- data.frame(Class = 1:nlevels(as.factor(DivType$Class)),
                                  ClassName = paste("Class", 1:nlevels(as.factor(DivType$Class))),
                                  ClassColor = 1 + 1:nlevels(as.factor(DivType$Class)),
                                  CompetenceColor = nlevels(as.factor(DivType$Class)) + 1 + 1:nlevels(as.factor(DivType$Class)))
          }
        }
      }
    }
      
  # --> edit 2014/09/11 #
  else if (div.parameter == "progeny")
    {
    cellsAtFirstImage <- find.Cells(InfoVertices$Images[ini], Cells, Divisions, Vertices, InfoVertices)
    oldest.ancestors <- vector() # used later
    RGBs <- list()
    VERTs <- list()
    color.first <- c(49, 112, 86)
    color.delta <- c(41, 62, 83)
    visual.dif <- 40
    for (cell in cellsAtFirstImage)
      {
      lineage <- which(cell == cellsAtFirstImage)
      vert <- as.numeric(na.omit(as.numeric(Cells[Cells$Cell == cell, -1])))
      VERTs[lineage] <- list(vert)
      if (lineage == 1)
        {
        RGBs[lineage] <- list(color.first)
        }
      else
        {
        completed.neighbors <- vector()
        for (completed.lineage in 1:(lineage-1))
          {
          if (length(which(vert %in% unlist(VERTs[completed.lineage]))) > 0) # at least one vertex is shared = cells are neighbors
            {
            completed.neighbors <- c(completed.neighbors, completed.lineage)
            }
          }
        if (length(completed.neighbors) == 0) # edit 2014/12/15
          {
          RGBs[lineage] <- list(color.first)
          }
        else
          {
          RGB <- color.first + color.delta * length(completed.neighbors)
          RGB <- RGB - 255*pmax(0, floor((RGB-1)/255))
          cleanColors <- T
          checkColorOfNeighbors <- F
          while (!checkColorOfNeighbors)
            {
            for (completed.neighbor in completed.neighbors)
              {
              differentR <- (RGB[1] < unlist(RGBs[completed.neighbor])[1] - visual.dif | RGB[1] > unlist(RGBs[completed.neighbor])[1] + visual.dif)
              differentG <- (RGB[2] < unlist(RGBs[completed.neighbor])[2] - visual.dif | RGB[2] > unlist(RGBs[completed.neighbor])[2] + visual.dif)
              differentB <- (RGB[3] < unlist(RGBs[completed.neighbor])[3] - visual.dif | RGB[3] > unlist(RGBs[completed.neighbor])[3] + visual.dif)
              if (!((differentR & differentG) | (differentR & differentB) | (differentG & differentB)))
                {
                RGB <- unlist(RGBs[completed.neighbor]) + rev(color.delta)
                RGB <- RGB - 255*pmax(0, floor((RGB-1)/255))
                cleanColors <- F
                break
                }
              else
                {
                cleanColors <- T
                }
              }
            if (cleanColors) { checkColorOfNeighbors <- T }
            }
          RGBs[lineage] <- list(RGB)
          }
        }
      }
    }
  # edit 2014/09/11 <-- #

  # Cells
  for (i in 1:nrow(Div))
    {
    # Cell ID
    cell <- as.character(Div$Cell[i])
    
    # Cell coordinates
    vert <- na.omit(as.numeric(Cells[Cells$Cell == cell, -1]))
    XY <- cell.coord(vert, VertX, VertY, Vertices)
    
    # Color
    if (div.parameter == "div&comp")
      {
      div.class <- Div$DivisionClass[i]
      comp.class <- Div$CompetenceClass[i]
      ifelse (is.na(div.class), color <- NA, ifelse (is.null(ClassName), color <- "olivedrab2", color <- ClassName$DivisionColor[ClassName$Class == div.class]))
      if (is.na(color)) # The color of division dominates over the color of competence
        {
        ifelse (is.na(comp.class), color <- NA, ifelse (is.null(ClassName), color <- "olivedrab", color <- ClassName$CompetenceColor[ClassName$Class == comp.class]))
        }
      }
    
    else if (div.parameter == "div")
      {
      div.class <- Div$DivisionClass[i]
      ifelse (is.na(div.class), color <- NA, ifelse (is.null(ClassName), color <- "olivedrab2", color <- ClassName$DivisionColor[ClassName$Class == div.class]))
      }
    
    else if (div.parameter == "comp")
      {
      comp.class <- Div$CompetenceClass[i]
      ifelse (is.na(comp.class), color <- NA, ifelse (is.null(ClassName), color <- "olivedrab", color <- ClassName$CompetenceColor[ClassName$Class == div.class]))
      }    

    # --> edit 2014/09/11 #
    else if (div.parameter == "progeny")
      {
      c(cellAtFirstImage, direct.mother) := find.OldestAncestor(cell, InfoVertices$Time[ini], Divisions)
      lineage <- which(cellAtFirstImage == cellsAtFirstImage)
      RGB <- unlist(RGBs[lineage])
      color <- rgb(RGB[1], RGB[2], RGB[3], max = 255)
      oldest.ancestors <- c(oldest.ancestors, cellAtFirstImage)
      }
    # edit 2014/09/11 <-- #

    else
      {
      # --> edit 2015/12/16 #
      ifelse (mini != 0 & mini == maxi & !is.na(d[Div$Cell == cell]),
              color <- Colors,
              color <- Colors[round((1 - mini*(length(Colors)-1)/(maxi-mini)) + (length(Colors)-1)/(maxi-mini) * d[Div$Cell == cell])])
      #color <- Colors[round((1 - mini*(length(Colors)-1)/(maxi-mini)) + (length(Colors)-1)/(maxi-mini) * d[Div$Cell == cell])]
      # edit 2015/12/16 <-- #
      }
   
    # Plot
    if (is.na(color)) { color <- col.font } # edit 2016/07/20
    polygon(XY, col = color, border = "grey")
    # --> edit 2016/01/15 #
    if (show.cell.number)
      {
      centro <- get.Centroid(XY)
      #text(centro, substr(cell, 6, nchar(cell)), cex = 0.5, col = col.font)
      text(centro, substr(cell, 6, nchar(cell)), cex = 0.5, col = ifelse (black, "black", "white")) # edit 2016/07/20
      }
    # edit 2016/01/15 <-- #
    }

  # --> edit 2014/09/11 and 2014/12/15 #
  if (div.parameter == "progeny")
    {
    for (oldest.ancestor in unique(oldest.ancestors))
      {
      vert <- na.omit(as.numeric(Cells[Cells$Cell == oldest.ancestor, -1]))
      ifelse (before, img <- Image1, img <- Image2)
      c(VerticesX, VerticesY) := process.Vertices(img, LeafShape, meanAngle, InfoVertices, Vertices)
      if (alignToPetioleLaminaBoundary)
        {
        Offset <- alignToPLB(VerticesY, cellAtPetioleLaminaBoundary, Cells, Vertices)
        VerticesY <- VerticesY - Offset
        }
      XY <- cell.coord(vert, VerticesX, VerticesY, Vertices)
      polygon(XY, col = NA, border = "grey", lwd = 2)
      }
    }
  # edit 2014/09/11 and 2014/12/15 <-- #

  # Scale bar
  if (!missing(scaleBar)) { plot.ScaleBar(ylim, scaleBar, tick, col.font) }
     
  # Scale of the division parameter (if required)
  if (growthScale)
    {
    if (div.parameter %in% c("CellArea", "AreaAtTime1", "AreaAtTime2", "Area_At_Division", "Symmetry", "Cell_Cycle_Duration"))
      {
      if (missing(txt.legend))
        {
        if (div.parameter == "CellArea" | (div.parameter == "AreaAtTime1" & before)) { txt.legend <- expression(paste("Cell area (µm"^2, ")")) }
        else if (div.parameter %in% c("AreaAtTime1", "AreaAtTime2")) { txt.legend <- expression(paste("Sector area (µm"^2, ")")) }
        else if (div.parameter == "Area_At_Division") { txt.legend <- expression(paste("Cell area at division (µm"^2, ")")) }
        else if (div.parameter == "Symmetry") { txt.legend <- "Symmetry of division" }
        else if (div.parameter == "Cell_Cycle_Duration") { txt.legend <- "Duration of cell cycle (h)" }
        }
      plot.GrowthScale(zlim, col.font, Percent = F, txt.legend = txt.legend, Colors = Colors, drawTicks = drawTicks, ticksGrowthScale = ticksGrowthScale)
      }
    #else
    else if (div.parameter %in% c("div&comp", "div", "comp")) # edit 2014/09/11
      {
      plot.DivAndCompLegend(div.parameter, col.font, ClassName)
      }
    # --> edit 2015/12/15 #
    else if (div.parameter == "Number_Of_Divisions")
      {
      if (missing(txt.legend)) { txt.legend <- "Number of divisions" }
      plot.DivisionScale(zlim, col.font, txt.legend, Colors)
      }
    # edit 2015/12/15 <-- #
    }
  }


plot.intervals <- function(Cells, Divisions, Vertices, InfoVertices, meanAngle,
                           alignToPetioleLaminaBoundary = T, cellAtPetioleLaminaBoundary,
                           interval, Image0, initial.stage = T, plot.type,
                           leafShape, leafShapeKriging = F, before = F, k.growth = 2,
                           PNG = T, black = T, xlim, ylim, scaleBar, tick = F,
                           Percent, round.zlim = T, zlim, fix.min = T, fix.max = T,
                           colorPaletteOfGFtbox = T, growthScale = T, drawTicks = T, txt.legend,
                           anisotropy, aniso.threshold = 0.05, aniso.lwd, aniso.lwd.constant = F,
                           n.pred = 100, polynomDegree = 3, rangeParameter = 0.001, contour.levels,
                           cellularScale = F, plotResidual = F, exportCellValues = F,
                           div.parameter = "div&comp", round.zlimDiv = T, zlimDiv, fix.minDiv = T, fix.maxDiv = T,
                           colorPaletteOfGFtboxDiv = T, growthScaleDiv = T, drawTicksDiv = F, txt.legendDiv,
                           wd, File, ini, show.cell.number = F)
  {
  # This function plots growth, kriging, bivariate interpolation or cell division at a time step 'interval'
  # given an initial 'Image 0'. It searches which images in the data set
  # are at the closest congruence to 'interval' until the end of the data set is reached.
  #
  # Arguments are essentially the same as for 'plot.growth()', 'plot.kriging()' and 'plot.division()', except that:
  #
  #   - growth is recomputed at each step, so 'Growth' and 'Shapes' are unuseful here.
  #
  #   - 'initial.stage' should be set to TRUE to draw:
  #       - cell outlines at 'Image0' if plot.type = "growth"
  #       - competence of all cells at 'Image0' if plot.type = "division"
  #
  #   - the time step 'interval' should be given in hours.
  #
  #   - 'plot.type' defines which of growth, kriging or division is to be plotted.
  #     One of this should be entered:
  #       - plot.type = "growth"
  #       - plot.type = "kriging"
  #       - plot.type = "division"
  #
  #   - 'PNG' is set to TRUE by default
  #
  #   - the colors for growth/kriging and division are independent, so the relevant arguments
  #     are duplicated with the suffix 'Div' for division
  #
  #   - 'wd' is the working directory
  #
  #   - 'File' is the name of the PointTracker raw csv file
  #
  #   - 'ini' is the index of the image from which divisions start to be processed (see 'process.AllDiv()').


  if (!plot.type %in% c("growth", "kriging", "division"))
    {
    stop("The argument 'plot.type' should be one string within \"growth\", \"kriging\" or \"division\"")
    }
  ifelse (plot.type == "kriging", plot.name <- "contours", plot.name <- plot.type)
  
  if (!file.exists("Graphical Outputs")) { dir.create("Graphical Outputs") }
  if (!file.exists("Divisions") & plot.type == "division") { dir.create("Divisions") } 
      
  if (plot.type %in% c("growth", "kriging"))
    {
    if      (k.growth == 2) { pref <- "karea" }
    else if (k.growth == 3) { pref <- "kmaj" }
    else if (k.growth == 4) { pref <- "kmin" }
    else if (k.growth == 5) { pref <- "theta" }
    else if (k.growth == 6) { pref <- "phi" }
    else if (k.growth == 7) { pref <- "anisotropy" }
    else if (k.growth == 8) { pref <- "kpertoml" }
    else if (k.growth == 9) { pref <- "kml" }
    else if (k.growth == 10) { pref <- "cellarea" }
    else if (k.growth == 11) { pref <- "sectorarea" }

    if (plot.type == "kriging")
      {
      if      (!cellularScale & !leafShapeKriging) { pref <- paste(pref, "TrackedRegion", sep = "__") }
      else if (!cellularScale & leafShapeKriging) { pref <- paste(pref, "WholeLeaf", sep = "__") }
      else if (cellularScale & !plotResidual) { pref <- paste(pref, "Cell", sep = "__") }
      else if (cellularScale & plotResidual) { pref <- paste(pref, "CellResidual", sep = "__") }
      }
    }
  else
    {
    if      (div.parameter == "div&comp") { pref <- "division and competence" }
    else if (div.parameter == "div") { pref <- "division only" }
    else if (div.parameter == "comp") { pref <- "competence only" }
    else if (div.parameter == "progeny") { pref <- "progeny" } # edit 2014/09/11
    else if (div.parameter == "CellArea") { pref <- "area of all cells" }
    else if (div.parameter %in% c("AreaTime1", "AreaTime2")) { pref <- "area of all sectors" }
    else if (div.parameter == "Area_At_Division") { pref <- "cell area at division" }
    else if (div.parameter == "Symmetry") { pref <- "symmetry of division" }
    else if (div.parameter == "Cell_Cycle_Duration") { pref <- "duration of cell cycle" } 
    else if (div.parameter == "Number_Of_Divisions") { pref <- "number of divisions" } # edit 2015/12/15
    }
      
  if (missing(Image0)) { Image0 <- as.character(InfoVertices$Images[1]) }
 
  time0 <- InfoVertices$Time[InfoVertices$Images == Image0]
  timeFinal <- InfoVertices$Time[nrow(InfoVertices)]
  wanted.intervals <- seq(time0, timeFinal, interval)
  if (wanted.intervals[length(wanted.intervals)] + interval < timeFinal + interval/3)
    {
    wanted.intervals <- c(wanted.intervals, wanted.intervals[length(wanted.intervals)] + interval)
    # check if the same image has not been incorporated twice
    actualTimeNminus1 <- InfoVertices$Time[which.min((InfoVertices$Time - wanted.intervals[length(wanted.intervals)-1])^2)]
    actualTimeN <- InfoVertices$Time[which.min((InfoVertices$Time - wanted.intervals[length(wanted.intervals)])^2)]
    ImageNminus1 <- as.character(InfoVertices$Image[InfoVertices$Time == actualTimeNminus1])
    ImageN <- as.character(InfoVertices$Image[InfoVertices$Time == actualTimeN])
    if (ImageNminus1 == ImageN) { wanted.intervals <- wanted.intervals[1:(length(wanted.intervals)-1)] }
    }
  if (initial.stage & plot.type %in% c("growth", "division")) { wanted.intervals <- c(wanted.intervals[1], wanted.intervals) }

  for (i in 1:(length(wanted.intervals)-1))
    {
    actualTime1 <- InfoVertices$Time[which.min((InfoVertices$Time - wanted.intervals[i])^2)]
    actualTime2 <- InfoVertices$Time[which.min((InfoVertices$Time - wanted.intervals[i+1])^2)]
    Image1 <- as.character(InfoVertices$Image[InfoVertices$Time == actualTime1])
    Image2 <- as.character(InfoVertices$Image[InfoVertices$Time == actualTime2])

    c(Growth, Shapes) := recompute.GrowthAndShapes(Cells, Divisions, Vertices, InfoVertices, Image1, Image2, meanAngle)
     
    if (PNG)
      {
      filename <- paste("Graphical Outputs/Intervals__", plot.name, "__", pref, "__",
                        round(wanted.intervals[i]), " - ", round(wanted.intervals[i+1]), " h__(",
                        substr(Image1, 1, nchar(Image1) - 4), " - ", substr(Image2, 1, nchar(Image2) - 4), ", ", 
                        round(actualTime1, 1), " - ", round(actualTime2, 1), " h).png", sep = "")
      if (nchar(filename) > 259 - nchar(wd)) { filename <- paste(substr(filename, 1, 259 - nchar(wd) - 4), ".png", sep = "") }
      size.in.pix <- 2000
      png(filename, width = size.in.pix, height = size.in.pix, res = 0.15*size.in.pix)
      }
    
    if (plot.type == "growth")
      {
      if (i == 1 & initial.stage)
        {
        plot.growth(InfoVertices, Vertices, Cells, Divisions, meanAngle, Growth, Shapes,
                    alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary,
                    Image1, Image2, leafShape, before, k.growth,
                    PNG, black, xlim, ylim, scaleBar, tick,
                    Percent, round.zlim, zlim, fix.min, fix.max,
                    colorPaletteOfGFtbox, growthScale = F, drawTicks, txt.legend,
                    anisotropy = F, aniso.threshold, aniso.lwd, aniso.lwd.constant)
        }
      else
        {
        plot.growth(InfoVertices, Vertices, Cells, Divisions, meanAngle, Growth, Shapes,
                    alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary,
                    Image1, Image2, leafShape, before, k.growth,
                    PNG, black, xlim, ylim, scaleBar, tick,
                    Percent, round.zlim, zlim, fix.min, fix.max,
                    colorPaletteOfGFtbox, growthScale, drawTicks, txt.legend,
                    anisotropy, aniso.threshold, aniso.lwd, aniso.lwd.constant)
        }
      }
    
    else if (plot.type == "kriging")
      {
      if (missing(anisotropy)) { anisotropy <- F }
      plot.kriging(InfoVertices, Vertices, Cells, Divisions, meanAngle, Growth, Shapes,
                   alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary,
                   Image1, Image2, leafShape, leafShapeKriging, before, k.growth,
                   PNG, black, xlim, ylim, scaleBar, tick,
                   Percent, round.zlim, zlim, fix.min, fix.max,
                   colorPaletteOfGFtbox, growthScale, drawTicks, txt.legend,
                   anisotropy, aniso.threshold, aniso.lwd, aniso.lwd.constant,
                   n.pred, polynomDegree, rangeParameter, contour.levels,
                   cellularScale, plotResidual, exportCellValues)
      }
    
    else if (plot.type == "division")
      {
      Div <- process.DivWithInterval(InfoVertices, Vertices, Cells, Divisions, meanAngle, Growth, Shapes,
                                     alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary,
                                     Image1, Image2, leafShape, before, File, ini)
      #if (i == 1 & initial.stage & div.parameter %in% c("Area_At_Division", "Symmetry", "Cell_Cycle_Duration"))
      if (i == 1 & initial.stage & div.parameter %in% c("Area_At_Division", "Symmetry", "Cell_Cycle_Duration", "Number_Of_Divisions")) # edit 2015/12/15
        {
        plot.division(InfoVertices, Vertices, Cells, Divisions, meanAngle, Shapes, Div, 
                      alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary,
                      Image1, Image2, leafShape, before,
                      PNG, black, xlim, ylim, scaleBar, tick,
                      div.parameter, round.zlimDiv, zlimDiv, fix.minDiv, fix.maxDiv,
                      colorPaletteOfGFtboxDiv, growthScale = F, drawTicksDiv, txt.legendDiv)
        }
      else
        {
        plot.division(InfoVertices, Vertices, Cells, Divisions, meanAngle, Shapes, Div, 
                      alignToPetioleLaminaBoundary, cellAtPetioleLaminaBoundary,
                      Image1, Image2, leafShape, before,
                      PNG, black, xlim, ylim, scaleBar, tick,
                      div.parameter, round.zlimDiv, zlimDiv, fix.minDiv, fix.maxDiv,
                      #colorPaletteOfGFtboxDiv, growthScaleDiv, drawTicksDiv, txt.legendDiv)
                      colorPaletteOfGFtboxDiv, growthScaleDiv, drawTicksDiv, txt.legendDiv, ini, show.cell.number) # edits 2014/09/11 and 2016/01/15
        }

      ifelse (i == 1, DivTotal <- Div, DivTotal <- rbind(DivTotal, Div))
      
      if (i == length(wanted.intervals)-1)
        {
        csvname <- paste("Divisions/Divisions__",
                        substr(Image0, 1, nchar(Image0) - 4), " - ", substr(Image2, 1, nchar(Image2) - 4), "__", 
                        interval, " h.csv", sep = "")
        write.csv(DivTotal, csvname, row.names = F)
        }
      
      }
    
    if (PNG) { dev.off() }
    }

  return (length(wanted.intervals)-1) # counts the plotted images
  }


draw.Lineage <- function (target.cell, ProcessedImages, InfoVertices, Divisions, Cells, Vertices,
                          leafShape, meanAngle, outputDir, n.col, n.row, main.plot, display.vertices,
                          makeAVI, compr, fps)
  {                       
  # This function plots the whole lineage (ancestors and descendants) of a target cell
  # and save the XY coordinates in 'outputDir'.
  # The target cell should be given as an integer.
  
  # Create or empty output directory
  if (!file.exists(outputDir)) { dir.create(outputDir) } 
  else if (length(dir(outputDir)) > 0) { file.remove(paste(outputDir, dir(outputDir), sep = "/")) }

  # Initial settings
  lineage <- find.Lineage(target.cell, Divisions)
  daughters <- cell <- paste("Cell", target.cell)
  level_mothers <- level_daughtersMax <- level_daughters <- 0
  
  # Find the maximal number of mothers
  i <- 1
  c(oldest.ancestor, direct.mother) := find.OldestAncestor(cell, InfoVertices$Time[i], Divisions)
  all.mothers <- temp.oldest.ancestor <- oldest.ancestor
  while (oldest.ancestor != cell)
    {
    i <- i + 1
    c(new.oldest.ancestor, direct.mother) := find.OldestAncestor(cell, InfoVertices$Time[i], Divisions)
    if (new.oldest.ancestor != oldest.ancestor)
      {
      oldest.ancestor <- new.oldest.ancestor
      level_mothers <- level_mothers + 1
      all.mothers <- c(all.mothers, oldest.ancestor)
      }
    }

  # Find the maximal number of daughters
  all.daughters <- vector()
  while (length(daughters) != 0)
    {
    for (d in daughters)
      {
      if (d %in% Divisions$Cell)
        {
        daughters <- c(daughters, paste("Cell", c(Divisions[Divisions$Cell == d, "Daughter cell 1"], Divisions[Divisions$Cell == d, "Daughter cell 2"])))
        level_daughtersMax <- level_daughtersMax + 1
        all.daughters <- c(all.daughters, paste("Cell", c(Divisions[Divisions$Cell == d, "Daughter cell 1"], Divisions[Divisions$Cell == d, "Daughter cell 2"])))
        }
      daughters <- daughters[-which(daughters == d)]
      }
    }
  daughters <- cell
    
  # Color settings
  visual.dif <- 30
  RGB_target <- list(c(255, 93, 239))
  RGB_mothers <- list(c(0, 150, 255))
  if (level_mothers > 1)
    {
    for (m in 2:level_mothers)
      {
      RGB_mothers[m] <- list(round(runif(3, min = 0, max = 255)))
      while (!((unlist(RGB_mothers[m])[1] < unlist(RGB_mothers[m-1])[1] - visual.dif | unlist(RGB_mothers[m])[1] > unlist(RGB_mothers[m-1])[1] + visual.dif) &
               (unlist(RGB_mothers[m])[2] < unlist(RGB_mothers[m-1])[2] - visual.dif | unlist(RGB_mothers[m])[2] > unlist(RGB_mothers[m-1])[2] + visual.dif) &
               (unlist(RGB_mothers[m])[3] < unlist(RGB_mothers[m-1])[3] - visual.dif | unlist(RGB_mothers[m])[3] > unlist(RGB_mothers[m-1])[3] + visual.dif)))
        {
        RGB_mothers[m] <- list(round(runif(3, min = 0, max = 255)))
        }
      #if (m %% 2 == 0) { RGB_mothers[m] <- list(c(0, 0, 255)) }
      #else { RGB_mothers[m] <- RGB_mothers[1] }
      }
    }
  RGB_daughters <- list(c(38, 139, 0))
  if (level_daughtersMax > 1)
    {
    for (d in 2:level_daughtersMax)
      {
      RGB_daughters[d] <- list(round(runif(3, min = 0, max = 255)))
      while (!((unlist(RGB_daughters[d])[1] < unlist(RGB_daughters[d-1])[1] - visual.dif | unlist(RGB_daughters[d])[1] > unlist(RGB_daughters[d-1])[1] + visual.dif) &
               (unlist(RGB_daughters[d])[2] < unlist(RGB_daughters[d-1])[2] - visual.dif | unlist(RGB_daughters[d])[2] > unlist(RGB_daughters[d-1])[2] + visual.dif) &
               (unlist(RGB_daughters[d])[3] < unlist(RGB_daughters[d-1])[3] - visual.dif | unlist(RGB_daughters[d])[3] > unlist(RGB_daughters[d-1])[3] + visual.dif)))
        {
        RGB_daughters[d] <- list(round(runif(3, min = 0, max = 255)))
        }
      #if (d %% 2 == 0) { RGB_daughters[d] <- list(c(93, 255, 109)) }
      #else { RGB_daughters[d] <- RGB_daughters[1] }
      }
    }

  # Loop on each image
  for (my_img in 1:nrow(ProcessedImages))
    {

    # Image properties
    img <- as.character(ProcessedImages[my_img, 1])
    idx_img <- which(InfoVertices$Images == img)
    Time <- InfoVertices$Time[idx_img]
    AngleShift <- -InfoVertices$AngleShift[idx_img]
    scaleX <- InfoVertices$ScalingX[idx_img]
    scaleY <- InfoVertices$ScalingY[idx_img]
    shiftX <- InfoVertices$ShiftX[idx_img]
    shiftY <- InfoVertices$ShiftY[idx_img]

    # Get leaf outline
    if (leafShape) { LeafShape <- get.leafShape(img, meanAngle, InfoVertices) }
    else { LeafShape <- data.frame(V1 = c(0, 0), V2 = c(0, 0), X = c(0, 0), Y = c(0, 0)) }
    
    # Check if the cell has divided at image 'img'
    isDaughter <- F
    if (cell %in% Divisions$Cell)
      {
      if (Divisions$Time[Divisions$Cell == cell] < idx_img)
        {
        isDaughter <- T
        # Check how many daughter cells have divided
        for (d in daughters)
          {
          if (d %in% Divisions$Cell)
            {
            if (Divisions$Time[Divisions$Cell == d] < idx_img)
              {
              daughters <- c(daughters, paste("Cell", c(Divisions[Divisions$Cell == d, "Daughter cell 1"], Divisions[Divisions$Cell == d, "Daughter cell 2"])))
              daughters <- daughters[-which(daughters == d)]
              level_daughters <- level_daughters + 1
              }
            }
          }
        }
      }

    # If the cell is not present at image 'img', find the oldest ancestor
    c(oldest.ancestor, direct.mother) := find.OldestAncestor(cell, Time, Divisions)
    if (cell == oldest.ancestor)
      {
      isMother <- F
      level_mothers <- 0
      }
    else
      {
      isMother <- T
      # Check how many generations apart is the mother
      if (oldest.ancestor != temp.oldest.ancestor)
        {
        temp.oldest.ancestor <- oldest.ancestor
        level_mothers <- level_mothers - 1
        }
      }
  
    # Coordinates
    vert <- na.omit(as.numeric(Cells[Cells$Cell == oldest.ancestor, -1]))
    VerticesX <- as.numeric(Vertices[, which(names(Vertices) == img)[1]])
    VerticesY <- as.numeric(Vertices[, which(names(Vertices) == img)[2]])
    VerticesX <- VerticesX / scaleX - shiftX
    VerticesY <- VerticesY / scaleY - shiftY
    rotated.vert <- fn.rotation(cbind(VerticesX, VerticesY), -AngleShift * pi/180 )
    VerticesX <- rotated.vert[,1]
    VerticesY <- rotated.vert[,2]
    XY <- cell.coord(vert, VerticesX, VerticesY, Vertices)
  
    # Color choice
    if (isMother) { c(R, G, B) := unlist(RGB_mothers[level_mothers]) }    
    else if (isDaughter) { c(R, G, B) := unlist(RGB_daughters[level_daughters]) }      
    else { c(R, G, B) := unlist(RGB_target[1]) }
      
    # Save information
    write.table(XY, paste(outputDir, "/cell_XY_", idx_img, ".txt", sep = ""), row.names = F, col.names = F, sep = "\t")
    write.table(data.frame(idx_img, AngleShift, scaleX, scaleY, shiftX, shiftY, isMother, isDaughter, img, paste(round(Time, 1), "h"),
                           level_mothers, level_daughters, R, G, B),
                paste(outputDir, "/Image_Info.txt", sep = ""), row.names = F, col.names = F, sep = "\t", append = T)
  
    if (my_img == 1)
      {
      cellMax <- oldest.ancestor
      write.table(data.frame(meanAngle, n.col, n.row, main.plot, display.vertices, makeAVI, compr, fps), paste(outputDir, "/Montage_Info.txt", sep = ""), row.names = F, col.names = F, sep = "\t")
      }
    
    if (my_img == nrow(ProcessedImages))
      {
      #vert <- na.omit(as.numeric(Cells[Cells$Cell == cellMax, -1]))
      vert <- na.omit(as.numeric(Cells[Cells$Cell == oldest.ancestor, -1]))
      VerticesX <- as.numeric(Vertices[, which(names(Vertices) == img)[1]])
      VerticesY <- as.numeric(Vertices[, which(names(Vertices) == img)[2]])
      VerticesX <- VerticesX / scaleX - shiftX
      VerticesY <- VerticesY / scaleY - shiftY
      rotated.vert <- fn.rotation(cbind(VerticesX, VerticesY), -AngleShift * pi/180 )
      VerticesX <- rotated.vert[,1]
      VerticesY <- rotated.vert[,2]
      XY <- cell.coord(vert, VerticesX, VerticesY, Vertices)
      write.table(XY, paste(outputDir, "/cell_XY_max.txt", sep = ""), row.names = F, col.names = F, sep = "\t")
      }

    # --> edit 2014/09/11 #
    if (my_img == 1)
      {
      plotCount <- 1
      if (leafShape) # Plot the cell/progeny within the leaf outline of the first image
        {
        LeafSpline <- get.leafSpline(LeafShape)
        xlim <- c(-1.01*max(abs(LeafSpline$x)), 1.01*max(abs(LeafSpline$x)))
        ylim <- c(min(LeafSpline$y) - 0.01*abs(min(LeafSpline$y)), 1.01*max(LeafSpline$y))
        plot(LeafShape$X, LeafShape$Y, asp = 1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = xlim, ylim = ylim)
        polygon(LeafSpline, lwd = 2, border = "grey")
        c(VerticesX, VerticesY) := process.Vertices(img, LeafShape, meanAngle, InfoVertices, Vertices)
        for (lin in c(all.mothers, all.daughters)[which(c(all.mothers, all.daughters) %in% find.Cells(img, Cells, Divisions, Vertices, InfoVertices))])
          {
          vert <- na.omit(as.numeric(Cells[Cells$Cell == lin, -1]))
          XY <- cell.coord(vert, VerticesX, VerticesY, Vertices)
          polygon(XY, border = "grey", lwd = 2)
          }
        }
      else # Empty plot
        {
        plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
        }
      }

    plotCount <- plotCount + 1
    if (n.row > 1 & plotCount %in% ((n.col+2)*(1:(n.row)-1))) # Empty plots
      {
      plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      plotCount <- plotCount + 2
      }
    # edit 2014/09/11 <-- #

    # Plot the cells
    #vert <- na.omit(as.numeric(Cells[Cells$Cell == cellMax, -1]))
    vert <- na.omit(as.numeric(Cells[Cells$Cell == oldest.ancestor, -1]))
    c(VerticesX, VerticesY) := process.Vertices(img, LeafShape, meanAngle, InfoVertices, Vertices)
    XY <- cell.coord(vert, VerticesX, VerticesY, Vertices)
    plot(XY, asp = 1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    if (main.plot == "NAME") { mtext(substitute(bold(imName), list(imName = substr(img, 1, nchar(img) - 4))), side = 1, cex = 0.75) }
    else if (main.plot == "NUMBER") { mtext(substitute(bold(imNumber), list(imNumber = idx_img)), side = 1, cex = 0.75) }
    else if (main.plot == "TIME") { mtext(substitute(bold(tm), list(tm = paste(round(Time, 1), "h"))), side = 1, cex = 0.75) }
  
    #for (lin in lineage[which(lineage %in% find.Cells(img, Cells, Divisions, Vertices, InfoVertices))])
    for (lin in c(all.mothers, all.daughters)[which(c(all.mothers, all.daughters) %in% find.Cells(img, Cells, Divisions, Vertices, InfoVertices))])
      {
      vert <- na.omit(as.numeric(Cells[Cells$Cell == lin, -1]))
      XY <- cell.coord(vert, VerticesX, VerticesY, Vertices)
      polygon(XY, border = "grey")
      centro <- get.Centroid(XY)
      text(centro, substr(lin, 6, nchar(lin)), cex = 0.75)
      }
    }

  # --> edit 2014/09/11 #
  if (leafShape) # Plot the cell/progeny within the leaf outline of the last image
    {
    plotCount <- plotCount + 1
    while (plotCount < (n.col+2)*n.row) # Empty plot(s), if any
      {
      plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      plotCount <- plotCount + 1
      }
    LeafSpline <- get.leafSpline(LeafShape)
    xlim <- c(-1.01*max(abs(LeafSpline$x)), 1.01*max(abs(LeafSpline$x)))
    ylim <- c(min(LeafSpline$y) - 0.01*abs(min(LeafSpline$y)), 1.01*max(LeafSpline$y))
    plot(LeafShape$X, LeafShape$Y, asp = 1, bty = "n", type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = xlim, ylim = ylim)
    polygon(LeafSpline, lwd = 2, border = "grey")
    c(VerticesX, VerticesY) := process.Vertices(img, LeafShape, meanAngle, InfoVertices, Vertices)
    for (lin in c(all.mothers, all.daughters)[which(c(all.mothers, all.daughters) %in% find.Cells(img, Cells, Divisions, Vertices, InfoVertices))])
      {
      vert <- na.omit(as.numeric(Cells[Cells$Cell == lin, -1]))
      XY <- cell.coord(vert, VerticesX, VerticesY, Vertices)
      polygon(XY, border = "grey", lwd = 2)
      }
    }
  # edit 2014/09/11 <-- #
  }
