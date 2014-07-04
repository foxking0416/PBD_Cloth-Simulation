// Cloth Simulation using Position Based Dynamics
// Copyright 2013 Xing Du

#include "constraint.h"
#include <cassert>

#ifndef EPSILON
#define EPSILON 0.00001f
#endif
using namespace glm;
//----------Constraint Class----------//
Constraint::Constraint() : 
    m_vertices(NULL),
    m_stiffness(1.0f)
{
   ;
}

Constraint::Constraint(ParticleList *verts, float stiff) : 
    m_vertices(verts),
    m_stiffness(stiff)
{
    ;
}

Constraint::Constraint(const Constraint& other) : 
    m_vertices(other.m_vertices),
    m_stiffness(other.m_stiffness)
{
    ;
}

Constraint::~Constraint()
{
    m_vertices = NULL;
}

bool Constraint::project_constraint()
{
    return true;
}

//----------FixedPointConstraint Class----------//
FixedPointConstraint::FixedPointConstraint() : 
    Constraint()
{
    ;
}

FixedPointConstraint::FixedPointConstraint(ParticleList *verts, unsigned int p0, const glm::vec3& fixedpoint) : 
    Constraint(verts, 1.0f),
    m_p0(p0),
    m_fixd_point(fixedpoint)
{
    ;
}

FixedPointConstraint::FixedPointConstraint(const FixedPointConstraint& other) : 
    Constraint(other),
    m_p0(other.m_p0),
    m_fixd_point(other.m_fixd_point)
{
    ;
}

FixedPointConstraint::~FixedPointConstraint()
{
    ;
}

bool FixedPointConstraint::project_constraint()
{// TODO: implement the project function for FixedPointConstraint.
    //return true if current position is OK. return false if the position is being projected.
    m_vertices->lock_pos(m_p0);
	glm::vec3 predictPos = m_vertices->predicted_pos(m_p0);
	
    //float value = 0.0f;
	float value = glm::length(predictPos - m_fixd_point);
    if(value < EPSILON)
        return true;

    glm::vec3 dp0 = m_fixd_point - predictPos;
    m_vertices->predicted_pos(m_p0) += dp0 * m_stiffness;
    return false;
}

//----------StretchConstraint Class----------//
StretchConstraint::StretchConstraint() : 
    Constraint()
{
    ;
}

StretchConstraint::StretchConstraint(ParticleList *verts, float stiff, unsigned int p1, unsigned int p2, float length) : 
    Constraint(verts, stiff),
    m_p1(p1),
    m_p2(p2),
    m_rest_length(length)
{
    ;
}

StretchConstraint::StretchConstraint(const StretchConstraint& other) : 
    Constraint(other),
    m_p1(other.m_p1),
    m_p2(other.m_p2),
    m_rest_length(other.m_rest_length)
{
    ;
}

StretchConstraint::~StretchConstraint()
{
    ;
}

bool StretchConstraint::project_constraint()
{// TODO: implement the project function for StretchConstraint.
    //return true if current position is OK. return false if the position is being projected.
    glm::vec3 p1, p2;
    p1 = m_vertices->predicted_pos(m_p1);
    p2 = m_vertices->predicted_pos(m_p2);

    float length = glm::length(p1 - p2);
    if(fabs(length - m_rest_length) < EPSILON)
        return true;

    glm::vec3 dp1, dp2;
	float w1 = m_vertices->inv_mass(m_p1);
	float w2 = m_vertices->inv_mass(m_p2);
	dp1 = -1 * w1 / (w1 + w2) * (length - m_rest_length) * (p1 - p2) / length;
	dp2 = w2 / (w1 + w2) * (length - m_rest_length) * (p1 - p2) / length;
    m_vertices->predicted_pos(m_p1) += dp1 * m_stiffness;
    m_vertices->predicted_pos(m_p2) += dp2 * m_stiffness;

    return false;
}

//----------BendConstraint Class----------//
BendConstraint::BendConstraint() : 
    Constraint()
{
    ;
}

BendConstraint::BendConstraint(ParticleList *verts, float stiff, unsigned int p1, unsigned int p2, unsigned int p3, unsigned int p4, float phi) : 
    Constraint(verts, stiff),
    m_p1(p1), m_p2(p2), m_p3(p3), m_p4(p4),
    m_phi(phi)
{
    ;
}

BendConstraint::BendConstraint(const BendConstraint& other) : 
    Constraint(other),
    m_p1(other.m_p1), m_p2(other.m_p2), m_p3(other.m_p3), m_p4(other.m_p4),
    m_phi(other.m_phi)
{
    ;
}

BendConstraint::~BendConstraint()
{
    ;
}

bool BendConstraint::project_constraint()
{// TODO: implement the project function for BendConstraint.
    //return true if current position is OK. return false if the position is being projected.
    glm::vec3 p1 = m_vertices->predicted_pos(m_p1),
              p2 = m_vertices->predicted_pos(m_p2),
              p3 = m_vertices->predicted_pos(m_p3),
              p4 = m_vertices->predicted_pos(m_p4);
	
	float w1 = m_vertices->inv_mass(m_p1);
	float w2 = m_vertices->inv_mass(m_p2);	
	float w3 = m_vertices->inv_mass(m_p3);
	float w4 = m_vertices->inv_mass(m_p4);

	vec3 newP1(0, 0, 0);
	vec3 newP2 = p2 - p1;
	vec3 newP3 = p3 - p1;	
	vec3 newP4 = p4 - p1;

	vec3 n1 = glm::cross(newP2, newP3) / length(glm::cross(newP2, newP3));
	vec3 n2 = glm::cross(newP2, newP4) / length(glm::cross(newP2, newP4));
	float d = dot(n1, n2);

	if(d > 1)
		d = 1;
	if(d < -1)
		d = -1;

	if(fabs(acos(d) - m_phi) <EPSILON)
		return true;

	vec3 q3 = (cross(newP2, n2) + cross(n1, newP2) * d) / length(cross(newP2, newP3));
	vec3 q4 = (cross(newP2, n1) + cross(n2, newP2) * d) / length(cross(newP2, newP4));
	vec3 q2 = -(cross(newP3, n2) + cross(n1, newP3) * d) / length(cross(newP2, newP3))
				-(cross(newP4, n1) + cross(n2, newP4) * d) / length(cross(newP2, newP4));
	vec3 q1 = -q2 - q3 - q4;

	float denominator = w1 * dot(q1,q1) + w2 * dot(q2,q2) + w3 * dot(q3,q3) + w4 * dot(q4,q4);

	vec3 dp1 = -1 * w1 * sqrt(1-d*d)*(acos(d) - m_phi) * q1 / denominator;
	vec3 dp2 = -1 * w2 * sqrt(1-d*d)*(acos(d) - m_phi) * q2 / denominator;
	vec3 dp3 = -1 * w3 * sqrt(1-d*d)*(acos(d) - m_phi) * q3 / denominator;
	vec3 dp4 = -1 * w4 * sqrt(1-d*d)*(acos(d) - m_phi) * q4 / denominator;

	/*
	glm::vec3 vecP1P2 = p2 - p1;
	glm::vec3 vecP1P3 = p3 - p1;
	glm::vec3 vecP1P4 = p4 - p1;
	glm::vec3 normal_p1p2p3 = glm::cross(vecP1P2, vecP1P3);
	glm::vec3 normal_p1p4p2 = glm::cross(vecP1P4, vecP1P2);
	float w3 = m_vertices->inv_mass(m_p3);
	float w4 = m_vertices->inv_mass(m_p4);

	float cosTheta = glm::dot(normal_p1p2p3, normal_p1p4p2) / glm::length(normal_p1p2p3) / glm::length(normal_p1p4p2);
	float angle = glm::acos(cosTheta);

	if(fabs(angle - m_phi) <EPSILON)
		return true;


	float t = ((p2[0] - p1[0]) * (p3[0] - p1[0]) + (p2[1] - p1[1]) * (p3[1] - p1[1]) + (p2[2] - p1[2]) * (p3[2] - p1[2])) / 
				(pow(p2[0] - p1[0], 2) + pow(p2[1] - p1[1], 2) + pow(p2[2] - p1[2], 2));

	vec3 midPt(p1[0] + t * (p2[0] - p1[0]), p1[1] + t * (p2[1] - p1[1]), p1[2] + t * (p2[2] - p1[2]));


	float deltaPhi3 = -w3 / (w3 + w4) * (angle - m_phi);
	float deltaPhi4 = w4 / (w3 + w4) * (angle - m_phi);

	vec3 vecP3Mid = p3 - midPt;
	vec3 vecP4Mid = p4 - midPt;


	float s3 = cos(deltaPhi3 / 2);
	float x3 = sin(deltaPhi3 / 2) * vecP1P2[0];
	float y3 = sin(deltaPhi3 / 2) * vecP1P2[1];
	float z3 = sin(deltaPhi3 / 2) * vecP1P2[2];
	mat3x3 rotationMatrix3(vec3(1 - 2*y3*y3 - 2*z3*z3, 2*x3*y3 + 2*s3*z3, 2*x3*z3 - 2*s3*y3),
							vec3(2*x3*y3 - 2*s3*z3, 1 - 2*x3*x3 - 2* z3*z3, 2*y3*z3 + 2*s3*x3),
							vec3(2*x3*z3 + 2*s3*y3, 2*y3*z3 - 2*s3*x3, 1-2*x3*x3 - 2*y3*y3));

	float s4 = cos(deltaPhi4 / 2);
	float x4 = sin(deltaPhi4 / 2) * vecP1P2[0];
	float y4 = sin(deltaPhi4 / 2) * vecP1P2[1];
	float z4 = sin(deltaPhi4 / 2) * vecP1P2[2];
	mat3x3 rotationMatrix4(vec3(1 - 2*y4*y4 - 2*z4*z4, 2*x4*y4 + 2*s4*z4, 2*x4*z4 - 2*s4*y4),
							vec3(2*x4*y4 - 2*s4*z4, 1 - 2*x4*x4 - 2* z4*z4, 2*y4*z4 + 2*s4*x4),
							vec3(2*x4*z4 + 2*s4*y4, 2*y4*z4 - 2*s4*x4, 1-2*x4*x4 - 2*y4*y4));


	vec3 newVecP3Mid = rotationMatrix3 * vecP3Mid;
	vec3 newVecP4Mid = rotationMatrix4 * vecP4Mid;
	vec3 newP3 = midPt + newVecP3Mid;
	vec3 newP4 = midPt + newVecP4Mid;*/

    //glm::vec3 dp1, dp2, dp3, dp4;

	//dp3 = newP3 - p3;
	//dp4 = newP4 - p4;


    m_vertices->predicted_pos(m_p1) += dp1 * m_stiffness;
    m_vertices->predicted_pos(m_p2) += dp2 * m_stiffness;
    m_vertices->predicted_pos(m_p3) += dp3 * m_stiffness;
    m_vertices->predicted_pos(m_p4) += dp4 * m_stiffness;

    return false;
}

//----------CollisionConstraint Class----------//
CollisionConstraint::CollisionConstraint() : 
    Constraint()
{
    ;
}

CollisionConstraint::CollisionConstraint(ParticleList *verts, unsigned int p0, const glm::vec3& q, const glm::vec3& n) : 
    Constraint(verts, 1.0f),
    m_p0(p0),
    m_ref_point(q),//intersection point
    m_normal(n)
{
    ;
}

CollisionConstraint::CollisionConstraint(const CollisionConstraint& other) : 
    Constraint(other),
    m_p0(other.m_p0),
    m_ref_point(other.m_ref_point),
    m_normal(other.m_normal)
{
    ;
}

CollisionConstraint::~CollisionConstraint()
{
    ;
}

bool CollisionConstraint::project_constraint()
{// TODO: implement the project function for CollisionConstraint.
    //return true if current position is OK. return false if the position is being projected.
    glm::vec3 p0 = m_vertices->predicted_pos(m_p0);
	glm::vec3 dir = p0 - m_ref_point;
	float value = glm::dot(dir, m_normal);
    //float value = 0.0f;
    if(value > 0.0f)
        return true;

    glm::vec3 dp0 = m_ref_point - p0;
    m_vertices->predicted_pos(m_p0) += dp0 * m_stiffness;

    return false;
}

//----------CollisionConstraint Class----------//
SelfCollisionConstraint::SelfCollisionConstraint() : 
    Constraint()
{
    ;
}

SelfCollisionConstraint::SelfCollisionConstraint(ParticleList *verts, unsigned int q, unsigned int p1, unsigned int p2, unsigned int p3, float h) :
    Constraint(verts, 1.0f),
    m_q(q), m_p1(p1), m_p2(p2), m_p3(p3),
    m_h(h)
{
    ;
}
SelfCollisionConstraint::SelfCollisionConstraint(const SelfCollisionConstraint& other) :
    Constraint(other),
    m_q(other.m_q), m_p1(other.m_p1), m_p2(other.m_p2), m_p3(other.m_p3),
    m_h(other.m_h)
{
    ;
}

SelfCollisionConstraint::~SelfCollisionConstraint()
{

}

bool SelfCollisionConstraint::project_constraint()
{
    glm::vec3 q, p1, p2, p3;
    q =  m_vertices->predicted_pos(m_q);
    p1 = m_vertices->predicted_pos(m_p1);
    p2 = m_vertices->predicted_pos(m_p2);
    p3 = m_vertices->predicted_pos(m_p3);

    q = q - p1;
    p2 = p2 - p1;
    p3 = p3 - p1;
    p1 = glm::vec3(0.0f);
    
    glm::vec3 normal(glm::cross(p2, p3));
    float c23 = glm::length(normal);
    normal = glm::normalize(normal);

    float value = glm::dot(q, normal) - m_h;
    if(value > 0.0f)
        return true;

    glm::vec3 dcq, dcp1, dcp2, dcp3;
    dcq = normal;
    dcp2 = (glm::cross(p3, q) + glm::cross(normal, p3) * glm::dot(normal, q)) / c23;
    dcp3 = -(glm::cross(p2, q) + glm::cross(normal, p2) * glm::dot(normal, q)) / c23;
    dcp1 = -dcq - dcp2 - dcp3;

    float wq, w1, w2, w3;
    wq = m_vertices->inv_mass(m_q);
    w1 = m_vertices->inv_mass(m_p1);
    w2 = m_vertices->inv_mass(m_p2);
    w3 = m_vertices->inv_mass(m_p3);

    float denominator = w1 * glm::dot(dcp1, dcp1) + w2 * glm::dot(dcp2, dcp2) + w3 * glm::dot(dcp3, dcp3) + wq * glm::dot(dcq, dcq);
    assert(denominator < EPSILON);

    glm::vec3 dq, dp1, dp2, dp3;
    float s = value / denominator;
    dq = -wq * s * dcq;
    dp1 = -w1 * s * dcp1;
    dp2 = -w2 * s * dcp2;
    dp3 = -w3 * s * dcp3;
    
    m_vertices->predicted_pos(m_q) += dq * m_stiffness;
    m_vertices->predicted_pos(m_p1) += dp1 * m_stiffness;
    m_vertices->predicted_pos(m_p2) += dp2 * m_stiffness;
    m_vertices->predicted_pos(m_p3) += dp3 * m_stiffness;
    return false;
}