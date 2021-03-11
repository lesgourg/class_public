Where to find information and documentation on CLASS?
======================================================

Author: Julien Lesgourgues


* __For what the code can actually compute__: all possible input parameters, all coded cosmological models, all functionalities, all observables, etc.: read the file `explanatory.ini` in the main `CLASS` directory: it is THE reference file where we keep track of all possible input and the definition of all input parameters. For that reason we recommend to leave it always unchanged and to work with copies of it, or with short input files written from scratch.


* __For the structure, style, and concrete aspects of the code__: this documentation, especially the `CLASS overview` chapter (the extensive automatically-generated part of this documentation is more for advanced users); plus the slides of our `CLASS` lectures, for instance those from New York 2019 available at

    `https://lesgourg.github.io/class-tour-NewYork.html`

    An updated overview of available `CLASS` lecture slides is always available at

    `http://lesgourg.github.io/courses.html`

    in the section `Courses on numerical tools`.


* __For the python wrapper of `CLASS`__: at the moment, the best are the "Usage I" and "Usage II" slides of the New York 2019 course,

    `https://lesgourg.github.io/class-tour-NewYork.html`

* __For the physics and equations used in the code__: mainly, the following papers:
    - *Cosmological perturbation theory in the synchronous and conformal Newtonian gauges*

     C. P. Ma and E. Bertschinger.

     http://arxiv.org/abs/astro-ph/9506072

     10.1086/176550

     Astrophys. J. __455__, 7 (1995)

    - *The Cosmic Linear Anisotropy Solving System (CLASS) II: Approximation schemes*

     D. Blas, J. Lesgourgues and T. Tram.

     http://arxiv.org/abs/1104.2933 [astro-ph.CO]

     10.1088/1475-7516/2011/07/034

     JCAP __1107__, 034 (2011)

    - *The Cosmic Linear Anisotropy Solving System (CLASS) IV: efficient implementation of non-cold relics*

     J. Lesgourgues and T. Tram.

     http://arxiv.org/abs/1104.2935 [astro-ph.CO]

     10.1088/1475-7516/2011/09/032

     JCAP __1109__, 032 (2011)

    - *Optimal polarisation equations in FLRW universes*

     T. Tram and J. Lesgourgues.

     http://arxiv.org/abs/1305.3261 [astro-ph.CO]

     10.1088/1475-7516/2013/10/002

     JCAP __1310__, 002 (2013)

    - *Fast and accurate CMB computations in non-flat FLRW universes*

     J. Lesgourgues and T. Tram.

     http://arxiv.org/abs/1312.2697 [astro-ph.CO]

     10.1088/1475-7516/2014/09/032

     JCAP __1409__, no. 09, 032 (2014)

    - *The CLASSgal code for Relativistic Cosmological Large Scale Structure*

     E. Di Dio, F. Montanari, J. Lesgourgues and R. Durrer.

     http://arxiv.org/abs/1307.1459 [astro-ph.CO]

     10.1088/1475-7516/2013/11/044

     JCAP __1311__, 044 (2013)

    - *The synergy between CMB spectral distortions and anisotropies*

     M. Lucca, N. Schöneberg, D. C. Hooper, J. Lesgourgues, J. Chluba.

     http://arxiv.org/abs/1910.04619 [astro-ph.CO]

     JCAP 02 (2020) 026

    - *Optimal Boltzmann hierarchies with nonvanishing spatial curvature*

     C. Pitrou, T. S. Pereira, J. Lesgourgues,

     http://arxiv.org/abs/2005.12119 [astro-ph.CO]

     Phys.Rev.D 102 (2020) 2, 023511

    plus also some latex notes on specific sectors:

    - *Equations for perturbed recombination*

     (can be turned on optionally by the user since v2.1.0)

     L. Voruz.

     http://lesgourg.github.io/class_public/perturbed_recombination.pdf

    - *PPF formalism in Newtonian and synchronous gauge*

     (used by default for the fluid perturbations since v2.6.0)

     T. Tram.

     http://lesgourg.github.io/class_public/PPF_formalism.pdf