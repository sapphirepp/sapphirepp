# Examples {#examples}

 1. @subpage parallel-shock  
    This example showcases the capabilities of @sapphire in the context of a
    parallel shock scenario. It provides a simple 1D simulation of a shock wave
    moving through a plasma, and introduces the use of runtime parameters in
    @sapphire.
 2. @subpage steady-state-parallel-shock  
    This example solves directly for the steady-state of the parallel shock scenario.
 3. @subpage gyro-advection  
    This example illustrates how an isotropic distribution of particles gyrates
    in a magnetic field, while being advected by a background plasma flow.
 4. @subpage closure  
    In this example the effect of truncating the expansion at $l_{\rm max}$ is explored.
 5. @subpage synchrotron-radiation  
    This example demonstrates the implementation of synchrotron radiation losses
    within the Vlasov–Fokker–Planck framework. It introduces the Landau–Lifshitz
    radiation–reaction term that accounts for energy loss of relativistic particles
    in a magnetic field. The setup is used to simulate the *steady-state parallel shock*
    example, enabling direct comparison between radiative and non-radiative cases.
 6. @subpage scattering-only  
    Serving as a comprehensive guide to @sapphire, this example is highly
    recommended for new developers. It focuses on using scattering to reach a
    solution with diminishing multipoles.
 7. @subpage convergence-study  
    In this advanced example, we derive an analytic solution for the system of
    equations solved in @sapphire in a special scenario. This is used to verify
    the accuracy of the numerical methods by performing a convergence study.

<div class="section_buttons">

| Previous                        |                              Next |
| :------------------------------ | --------------------------------: |
| [Visualization](#visualization) | [Implementation](#implementation) |

</div>
