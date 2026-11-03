source('R/00_load_all.R');bd_load_all()
L <- matrix(0L,70,70);L[33:38,32:40] <- 19L; L[9:10,8:11] <- 19L
before <- .bd_hash(L)
a <- bd_sky_exclusion(L,0L)
b <- bd_sky_exclusion(L,4L)
stopifnot(identical(a,L>0),all(b[a]),sum(b)>sum(a),
          isTRUE(b[33,28]), !b[33,26],identical(before,.bd_hash(L)))
# Inject faint wings into a sky-only field; dilation should improve sky baseline.
x <- outer(seq_len(70),rep(1L,70));y <- outer(rep(1L,70),seq_len(70))
wing <- 3*exp(-((x-35.5)^2+(y-36)^2)/120)
noise <- .09*sin(x*.19+y*.43) + .06*cos(x*.33-y*.17)
image <- wing + noise
s0 <- bd_sky_noise(image,L,687,sky_mask_radius_px=0L)
s4 <- bd_sky_noise(image,L,687,sky_mask_radius_px=4L)
stopifnot(s4$background_n<s0$background_n,
    identical(s4$original_L_sha256,before),
    abs(s4$residual_sky)<=abs(s0$residual_sky)+.03)
cat('SKY_ONLY_DILATION_PASS\n')
