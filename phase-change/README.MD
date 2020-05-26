# Phase change models in Basilisk solver #

## Sharp interface models
This phase change model includes:  
- the mass flux calcualtion  
    - [reference](https://github.com/GabrielGLK/thesis-cases/blob/master/reference-papers/phase-change/mass-transfer-model/sharp-interface-model/simpified-sharp-model/Sun%20%E7%AD%89%E3%80%82%20-%202012%20-%20Development%20of%20a%20vapor%E2%80%93liquid%20phase%20change%20model%20f.pdf)
    - [code](https://github.com/GabrielGLK/phase-change/blob/master/phase-change/phase_change_models/sharp-interface-model/film-boiling/phase_change_code/phase-change.h#L143)
    - Feature:  
The method to calcualte temperature gradient is different, we use similar method by [Quentin](http://basilisk.fr/sandbox/qmagdelaine/phase_change/elementary_body.h).

- the vof equation
    - reconstruct one extended divergence-free entire domain
        - [reference](https://github.com/GabrielGLK/thesis-cases/blob/master/reference-papers/phase-change/VOF/extended-VOF/Malan%20%E7%AD%89%E3%80%82%20-%202020%20-%20A%20geometric%20vof%20method%20for%20interface%20resolved%20phas.pdf)
        - [code](https://github.com/GabrielGLK/phase-change/blob/master/phase-change/phase_change_models/sharp-interface-model/film-boiling/small_case.c#L223)
    - vof advection
        - [Leon Malan method](https://github.com/GabrielGLK/thesis-cases/blob/master/reference-papers/phase-change/VOF/extended-VOF/Malan%20%E7%AD%89%E3%80%82%20-%202020%20-%20A%20geometric%20vof%20method%20for%20interface%20resolved%20phas.pdf), [code](https://github.com/GabrielGLK/phase-change/blob/master/phase-change/phase_change_models/sharp-interface-model/film-boiling/small_case.c#L146)
        - [normal method](https://github.com/GabrielGLK/phase-change/blob/master/phase-change/phase_change_models/sharp-interface-model/film-boiling/small_case.c#L153)

- the energy equation
    - ghost cell method
        - [reference](https://github.com/GabrielGLK/thesis-cases/blob/master/reference-papers/phase-change/mass-transfer-model/sharp-interface-model/Zhang%20%E5%92%8C%20Ni%20-%202018%20-%20Direct%20numerical%20simulations%20of%20incompressible%20mul.pdf)
        - [code]
        (https://github.com/GabrielGLK/phase-change/blob/master/phase-change/phase_change_models/sharp-interface-model/film-boiling/small_case.c#L170)
    - [tracer.h](http://basilisk.fr/src/tracer.h) + [diffusion.h](http://basilisk.fr/src/diffusion.h)

- Results
    - Film boiling
