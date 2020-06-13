mgstats no_flux_diffusion (scalar tr, scalar f, double dt) {

  /**
  We allocate fields for the *volume correction*, the face fraction and the
  weighted diffusion coefficient. */

  scalar volume_correction[];
  face vector f_f[], diffusion_coefficient[];

  #if TREE && !AXI
    volume_correction.prolongation = volume_correction.refine = fraction_refine;
  #endif

  #if TREE && !AXI && FACE_FRACTION_REFINE // seems better without
    foreach_dimension() {
      diffusion_coefficient.x.prolongation = fraction_refine;
      diffusion_coefficient.x.refine = fraction_refine;
    }
  #endif
  
  /**
  If the tracer is associatied to the phase $f=0$ instead of the phase $f=1$,
  the *volume correction* becomes $1 - f$ and the face fraction has to be
  replaced by $1 - f_f$. */ 

  foreach() {
    f[] = clamp(f[], 0., 1.);
    volume_correction[] = cm[]*max(tr.inverse ? 1. - f[] : f[], F_ERR);
  }
  boundary ({f, tr, volume_correction});

  /**
  To compute the face fraction, we use the function of Jose-Maria Lopez Herrera
  defined in [fracface.h](/sandbox/lopez/fracface.h). */

  face_fraction (f, f_f);
  foreach_face()
    diffusion_coefficient.x[] = tr.D*fm.x[]*(tr.inverse ? 1. - f_f.x[] : f_f.x[]);
  boundary((scalar *){diffusion_coefficient});

  /**
  The diffusion equation is solved thanks to [diffusion.h](/src/diffusion.h): */

  return diffusion (tr, dt, D = diffusion_coefficient, theta = volume_correction);
}