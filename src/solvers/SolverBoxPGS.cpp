#include "solvers/SolverBoxPGS.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

#include <Eigen/Dense>


SolverBoxPGS::SolverBoxPGS(RigidBodySystem* _rigidBodySystem) : Solver(_rigidBodySystem)
{

}

// Use the exact signature from your version 1
void accumulateCoupledContacts(const Contact* contact, const JBlock& JMMinv, const RigidBody* body, Eigen::Vector3f& x) {
    const auto& neighbour_contacts = body->contacts;
    for (const auto& nb_ct : neighbour_contacts) {
        if (nb_ct == contact) continue;
        const JBlock& J_other = (nb_ct->body0 == body ? nb_ct->J0 : nb_ct->J1);
        x -= JMMinv * J_other.transpose() * nb_ct->lambda;
    }
}

void SolverBoxPGS::solve(float h)
{
    // TODO Implement a PGS method that solves for the
    //      contact constraint forces in @a rigidBodySystem. 
    //      Assume a boxed LCP formulation.
    //
    float hinv = 1.0f / h;
    std::vector<Contact*>& contacts = m_rigidBodySystem->getContacts();
    const int numContacts = contacts.size();

    // Build array of 3x3 diagonal matrices, one for each contact.
    // 
    std::vector<Eigen::Matrix3f> Acontactii;
    if( numContacts > 0 )
    {
        // Pre-allocate and pre-compute to align with your requirement
        Acontactii.resize(numContacts);
        std::vector<Eigen::Vector3f> b(numContacts);

        // TODO Compute the right-hand side vector : b = -gamma*phi/h - J*vel - dt*JMinvJT*force
        //
        for (int i = 0; i < numContacts; i++) {
            auto& ct = contacts[i];
            
            // Using direct member access (xdot, omega, f, tau) as per your RigidBody class
            b[i] =
                -(ct->J0.block<3, 3>(0, 0) * ct->body0->xdot + ct->J0.block<3, 3>(0, 3) * ct->body0->omega)
                - (ct->J1.block<3, 3>(0, 0) * ct->body1->xdot + ct->J1.block<3, 3>(0, 3) * ct->body1->omega)
                - h * (ct->J0Minv.block<3, 3>(0, 0) * ct->body0->f + ct->J0Minv.block<3, 3>(0, 3) * ct->body0->tau)
                - h * (ct->J1Minv.block<3, 3>(0, 0) * ct->body1->f + ct->J1Minv.block<3, 3>(0, 3) * ct->body1->tau);
            
            float gamma = h * ct->k / (h * ct->k + ct->b);
            b[i] -= gamma * ct->phi * hinv;

            // TODO Compute the diagonal term : Aii = J0*Minv0*J0^T + J1*Minv1*J1^T
            // Pre-computing Aii outside the iteration loop for efficiency
            Acontactii[i] = ct->J0Minv * ct->J0.transpose() + ct->J1Minv * ct->J1.transpose();
            
            // Add regularization / CFM
            Acontactii[i](0, 0) += 1.0f / (ct->k * h * h + h * ct->b);
            // Optionally add to tangents if needed, but keeping it to (0,0) as per your version 1
            
            // Initialize lambda
            ct->lambda.setZero();
        }


        // PGS main loop.
        // There is no convergence test here. Simply stop after @a maxIter iterations.
        //
        for(int iter = 0; iter < m_maxIter; ++iter)
        {
            // TODO For each contact, compute an updated value of contacts[i]->lambda
            //      using matrix-free pseudo-code provided in the course appendix.
            //
            for(int i = 0; i < numContacts; ++i)
            {
                // TODO initialize current solution as x = b(i)
                //
                Contact* ct = contacts[i];
                Eigen::Vector3f x = b[i];


                // TODO Loop over all other contacts involving c->body0
                //      and accumulate : x -= (J0*Minv0*Jother^T) * lambda_other
                //
                accumulateCoupledContacts(ct, ct->J0Minv, ct->body0, x);


                // TODO Loop over all other contacts involving c->body1
                //      and accumulate : x -= (J0*Minv0*Jother^T) * lambda_other
                //
                accumulateCoupledContacts(ct, ct->J1Minv, ct->body1, x);


                // TODO Update lambda by solving the sub-problem : Aii * contacts[i].lambda = x
                //      For non-interpenetration constraints, lambda is non-negative and should be clamped
                //      to the range (0... infinity).
                //      Otherwise friction constraints should be clamped to the range (-mu*lambda_n, mu*lambda_n)
                //      where lambda_n is the impulse of the corresponding non-interpentration constraint.
                //
                // Using the pre-computed Acontactii[i]
                ct->lambda = Acontactii[i].inverse() * x;
                
                if (ct->lambda[0] < 0) ct->lambda[0] = 0;
                float fric_clamp = ct->mu * ct->lambda[0];
                ct->lambda[1] = std::max(std::min(ct->lambda[1], fric_clamp), -fric_clamp);
                ct->lambda[2] = std::max(std::min(ct->lambda[2], fric_clamp), -fric_clamp);

            }
        }
    }
}