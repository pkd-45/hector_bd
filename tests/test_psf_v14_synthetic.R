# Requires all dependencies of recovered v13 (Rfits, ProFound, ProFit, magicaxis).
source('R/psf_v14.R')
x<-c(88,89,-89,-88,20)
r<-.v14_axial_pa_clip(x)
stopifnot(all(r[1:4]),!r[5])
stopifnot(abs(stats::mad(c(-1,0,1))-1.4826)<1e-4)
for (d in c(21L,31L,61L)) {
  ee <- .v14_moffat_circle_ee(3,2.1,0.8,d)
  stopifnot(is.finite(ee),ee>0,ee<1)
}
stopifnot(.v14_moffat_circle_ee(3,2.1,0.8,41) >
          .v14_moffat_circle_ee(3,2.1,0.8,21))
p<-.build_default_psf(2.8,2.1,psf_dim=c(75L,75L))
stopifnot(nrow(p)>=21L,identical(attr(p,'ee_method'),
  'ELLIPTICAL_MOFFAT_ANALYTIC_CIRCULAR_APERTURE_720_ANGLES_APPROX'))
cat('D5 PSF synthetic helpers PASS\n')
