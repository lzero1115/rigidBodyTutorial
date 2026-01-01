#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"

Contact::Contact() : p(), n(), t1(), t2(), mu(0.4f),
    body0(nullptr), body1(nullptr), k(1e6f), b(1e5f), index(-1), pene(0.0f)
{

}

Contact::Contact(RigidBody* _body0, RigidBody* _body1, const Eigen::Vector3f& _p, const Eigen::Vector3f& _n, float _pene) :
    p(_p), n(_n), t1(), t2(), mu(0.4f), pene(_pene), body0(_body0), body1(_body1), k(1e6f), b(1e5f), index(-1)
{
    J0.setZero(3, 6);
    J1.setZero(3, 6);
    J0Minv.setZero(3, 6);
    J1Minv.setZero(3, 6);
    lambda.setZero(3);
    phi.setZero(3);
    phi(0) = _pene;

    body0->contacts.push_back(this);
    body1->contacts.push_back(this);
}

Contact::~Contact()
{

}

void Contact::computeContactFrame()
{
    // Compute the contact frame, which consists of an orthonormal
    //  bases formed the vector n, t1, and t2
    // 
    //  The first bases direction is given by the normal, n.
    //  Use it to compute the other two directions.

    // TODO Compute first tangent direction t1
    
    Eigen::Vector3f aux(1, 0, 0);
    if (std::abs(n.x()) > 0.9f) {
        aux = Eigen::Vector3f(0, 1, 0);
    }
    t1 = n.cross(aux).normalized();
    t2 = n.cross(t1).normalized();

}

void Contact::computeJacobian()
{
    // TODO Compute the Jacobians J0 and J1 
    // for body0 and body1, respectively.
    Eigen::Vector3f r0 = p - body0->x; 
    Eigen::Vector3f r1 = p - body1->x;

    J0.block<1, 3>(0, 0) = n.transpose();  
    J0.block<1, 3>(1, 0) = t1.transpose(); 
    J0.block<1, 3>(2, 0) = t2.transpose(); 
    // (w x r)*n = (r x n)*w
    J0.block<1, 3>(0, 3) = r0.cross(n).transpose();  
    J0.block<1, 3>(1, 3) = r0.cross(t1).transpose(); 
    J0.block<1, 3>(2, 3) = r0.cross(t2).transpose(); 

    J1.block<1, 3>(0, 0) = -n.transpose();
    J1.block<1, 3>(1, 0) = -t1.transpose();
    J1.block<1, 3>(2, 0) = -t2.transpose();

    J1.block<1, 3>(0, 3) = -(r1.cross(n).transpose());
    J1.block<1, 3>(1, 3) = -(r1.cross(t1).transpose());
    J1.block<1, 3>(2, 3) = -(r1.cross(t2).transpose());

    // Compute the J M^-1 blocks for each body. The code is provided.
    // 
    // However, together with the contact Jacobians J0 and J1, these will
    //   be used by the solver to assemble the blocked LCP matrices.
    //
    J0Minv.block(0,0,3,3) = (1.0f/body0->mass) * J0.block(0, 0, 3, 3);
    J0Minv.block(0,3,3,3) = J0.block(0, 3, 3, 3) * body0->Iinv;
    J1Minv.block(0,0,3,3) = (1.0f/body1->mass) * J1.block(0, 0, 3, 3);
    J1Minv.block(0,3,3,3) = J1.block(0, 3, 3, 3) * body1->Iinv;
}
