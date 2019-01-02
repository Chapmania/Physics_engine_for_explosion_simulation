using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public struct Collision
{
    public RigidBodyClass obj1, obj2;

    public Vector3 CollisionLocation;
    public Vector3 NormalDirection;
    public Vector3 Edge1Direction;
    public Vector3 Edge2Direction;
     
    public bool IsVertexFaceContact;

    public Collision(RigidBodyClass a1, RigidBodyClass b1, Vector3 p1, Vector3 n1, Vector3 ea1, Vector3 eb1, bool vf1)
    {
        obj1 = a1;
        obj2 = b1;
        CollisionLocation = p1;
        NormalDirection = n1;
        Edge1Direction = ea1;
        Edge2Direction = eb1;
        IsVertexFaceContact = vf1;
    }

}
