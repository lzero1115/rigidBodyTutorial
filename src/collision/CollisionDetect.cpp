#include "collision/CollisionDetect.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

CollisionDetect::CollisionDetect(RigidBodySystem* rigidBodySystem) : m_rigidBodySystem(rigidBodySystem)
{

}

void CollisionDetect::detectCollisions()
{
    // First, clear any existing contacts.
    //
    clear();

    // Next, loop over all pairs of bodies and test for contacts.
    //
    auto bodies = m_rigidBodySystem->getBodies();
    for(unsigned int i = 0; i < bodies.size(); ++i)
    {
        for(unsigned int j = i+1; j < bodies.size(); ++j)
        {
            RigidBody* body0 = bodies[i];
            RigidBody* body1 = bodies[j];

            // Special case: skip tests for pairs of static bodies.
            //
            if (body0->fixed && body1->fixed) 
                continue;

            // Test for sphere-sphere collision.
            if( body0->geometry->getType() == kSphere &&
                body1->geometry->getType() == kSphere )
            {
                collisionDetectSphereSphere(body0, body1);
            }
            // Test for sphere-box collision
            else if( body0->geometry->getType() == kSphere &&
                     body1->geometry->getType() == kBox )
            {
                collisionDetectSphereBox(body0, body1);
            }
            // Test for box-sphere collision (order swap)
            else if( body1->geometry->getType() == kSphere &&
                     body0->geometry->getType() == kBox )
            {
                collisionDetectSphereBox(body1, body0);
            }
        }
    }
}

void CollisionDetect::computeContactJacobians()
{
    // Build constraint Jacobians for all contacts
    //
    for(auto c : m_contacts)
    {
        c->computeContactFrame();
        c->computeJacobian();
    }
}

void CollisionDetect::clear()
{
    // First, remove all contacts from rigid bodies.
    //
    auto bodies = m_rigidBodySystem->getBodies();
    for (auto b : bodies)
    {
        b->contacts.clear();
    }

    // Then, cleanup the local contact array.
    //
    for(auto c : m_contacts)
    {
        delete c;
    }
    m_contacts.clear();

}

void CollisionDetect::collisionDetectSphereSphere(RigidBody* body0, RigidBody* body1)
{
    Sphere* sphere0 = dynamic_cast<Sphere*>(body0->geometry.get());
    Sphere* sphere1 = dynamic_cast<Sphere*>(body1->geometry.get());

    // Implement sphere-sphere collision detection.
    // The function should check if a collision exists, and if it does
    // compute the contact normal, contact point, and penetration depth.
    //
    Eigen::Vector3f vec = body0->x - body1->x;

    const float rsum = (sphere0->radius + sphere1->radius);
    const float dist = vec.norm();
    if( dist < rsum )
    {
        const Eigen::Vector3f n = vec / dist;
        const Eigen::Vector3f p = 0.5f * ((body0->x - sphere0->radius*n) + (body1->x + sphere1->radius*n));
        const float phi = dist-rsum;

        m_contacts.push_back( new Contact(body0, body1, p, n, phi) );
    }
}

void CollisionDetect::collisionDetectSphereBox(RigidBody* body0, RigidBody* body1)
{
    // TODO Implement sphere-box collision detection.
    //      The function should check if a collision exists.
    // 
    //      If it does, compute the contact normal, contact point, and penetration depth and
    //      create a Contact and add it to m_contacts.
    //

    Sphere* sphere = dynamic_cast<Sphere*>(body0->geometry.get());
    Box* box = dynamic_cast<Box*>(body1->geometry.get());

    Eigen::Vector3f c_local = body1->q.inverse() * (body0->x - body1->x);
    Eigen::Vector3f h = box->dim * 0.5f;

    // lazy clamp
    Eigen::Vector3f g;
    g.x() = std::max(-h.x(), std::min(h.x(), c_local.x()));
    g.y() = std::max(-h.y(), std::min(h.y(), c_local.y()));
    g.z() = std::max(-h.z(), std::min(h.z(), c_local.z()));

    Eigen::Vector3f delta = g - c_local;
    float distSq = delta.squaredNorm();
    float radius = sphere->radius;

    if (distSq < radius * radius){
        float dist = std::sqrt(distSq);
        Eigen::Vector3f n_local;
        float phi;

        if (dist > 1e-6f)
        {
            phi = dist - radius;
            n_local = - delta / dist; 
        }
        else // deep penetration, sphere center is inside box
        {
            float min_dist = std::numeric_limits<float>::max();
            int best_axis = 0;
            float side = 1.0f;

            for (int i = 0; i < 3; ++i)
            {
                // Distance to positive face and negative face
                float d_pos = h(i) - c_local(i);
                float d_neg = c_local(i) + h(i);

                if (d_pos < min_dist) { min_dist = d_pos; best_axis = i; side = 1.0f; }
                if (d_neg < min_dist) { min_dist = d_neg; best_axis = i; side = -1.0f; }
            }

            phi = -min_dist - radius; // Depth is distance to face + radius
            n_local.setZero();
            n_local(best_axis) = side;
            
            // Adjust 'g' to be on the closest face
            g = c_local;
            g(best_axis) = side * h(best_axis);
    
        }
        
        Eigen::Vector3f n_world = body1->q * n_local;
        Eigen::Vector3f p_world = body1->q * g + body1->x;
        
        m_contacts.push_back(new Contact(body0, body1, p_world, n_world, phi));

    }
}
