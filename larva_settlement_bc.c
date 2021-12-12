/* Trap boundary condition based on particle shear stress */

/* reflection component based on 2.5.1.3 example 1 in theory guide, 
which is the standard fluent boundary condition without an allowance for moving walls */


/* f: floor
 * dim: 2
 */

#include "udf.h"
#include  "dpm_mem.h"
#include "mem.h"

DEFINE_DPM_BC(larva_settle,tp,t,f,f_normal,dim)
{
    /* Variable declarations*/
    real alpha; /* angle of particle path with face normal */
    real vn = 0.;
    real nor_coeff = 1.; /* ideal reflection for normal velocity coefficient */
    real tan_coeff = 0.3; /* tanjential component is damped */
    real normal[3];
    real tau; /* shear stress*/
    cell_t c;
    int i, idim = dim;
    /*int zone_ID = THREAD_ID(t)*/
    real NV_VEC(x);

    #if RP_2D
    /* dim is always 2 in 2D compilation. Need special treatment for 2d
    axisymmetric and swirl flows */

    if (rp_axi_swirl)
    {
        real R = sqrt(TP_POS(tp)[1]*TP_POS(tp)[1] +
        TP_POS(tp)[2]*TP_POS(tp)[2]);
	    if (R > 1.e-20)
    	{
		idim = 3;
		normal[0] = f_normal[0];
		normal[1] = (f_normal[1]*TP_POS(tp)[1])/R;
		normal[2] = (f_normal[1]*TP_POS(tp)[2])/R;
	    }
	    else
        {
	        for (i=0; i<idim; i++)
                normal[i] = f_normal[i];
	    }
    }
    else
    #endif

    for (i=0; i<idim; i++)
        normal[i] = f_normal[i];


    if(TP_TYPE(tp) == DPM_TYPE_INERT) /* Only for inert particles */
    {
        /* get cell where the tracked particle is */
        tau = C_STORAGE_R(c,t,SV_WALL_SHEAR); 

        if(tau > 0.1) /* REFLECT IF SHEAR STRESS HIGHER THAN THE CRITICAL SHEAR STRESS - here 0.1Pa */
        {
            /* Calculate angle of incidence */
	        alpha = M_PI/2. - acos(MAX(-1.,MIN(1.,NV_DOT(normal,TP_VEL(tp))/
		        MAX(NV_MAG(TP_VEL(tp)),DPM_SMALL))));
	        if ((NNULLP(t)) && (THREAD_TYPE(t) == THREAD_F_WALL))
		        F_CENTROID(x,f,t);

	        /* calculate the normal component, rescale its magnitude by
	        the coefficient of restitution and subtract the change */

	        /* Compute normal velocity. */
	        for(i=0; i<idim; i++)
	            vn += TP_VEL(tp)[i]*normal[i];

	        /* Subtract off normal velocity. */
	        for(i=0; i<idim; i++)
	            TP_VEL(tp)[i] -= vn*normal[i];

	        /* Apply tangential coefficient of restitution. */
	        for(i=0; i<idim; i++)
	            TP_VEL(tp)[i] *= tan_coeff;

	        /* Add reflected normal velocity. */
	        for(i=0; i<idim; i++)
	            TP_VEL(tp)[i] -= nor_coeff*vn*normal[i];

	        /* Store new velocity in TP_VEL0 of particle */
	        for(i=0; i<idim; i++)
	            TP_VEL0(tp)[i] = TP_VEL(tp)[i];

	        return PATH_ACTIVE;
        }
        else
        {
            return PATH_ABORT; /* TRAP IF LESS THAN CRITICAL SHEAR STRESS */
        }
    }

return PATH_ABORT;

}





