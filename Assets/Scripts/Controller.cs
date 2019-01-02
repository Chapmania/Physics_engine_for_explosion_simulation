using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Controller : MonoBehaviour {

    List<RigidBodyClass> rigidBodies = new List<RigidBodyClass>();
    List<GameObject> blockList = new List<GameObject>();


    List<Collision> collisions = new List<Collision>();             //All collisions
    List<Collision> restingContact = new List<Collision>();         //Collisions involving Resting Contact
    List<Collision> collidingContact = new List<Collision>();       //Collisions involving Colliding Contact

    [SerializeField]
    GameObject TestMarker;
    bool test = false;

      // Use this for initialization
    void Start () {
		collisions.Clear();
		collidingContact.Clear();
		restingContact.Clear();
		CreateRigidBodyList();
	}
	
	// Update is called once per frame
	void Update () {
		collisions.Clear();
		collidingContact.Clear();
		restingContact.Clear();
		DetectCollisions();
        DefineCollisionType();
        CollidingContact();
	}

    void CreateRigidBodyList()
    {
        //Assumption : blockList[i] has same rigidBody as rigidBodies[i]. 
        GameObject[] blockArray = GameObject.FindGameObjectsWithTag("Block");
        for (int i = 0; i < blockArray.Length; i++)
        {
            blockList.Add(blockArray[i]);
            rigidBodies.Add(blockList[i].GetComponent<BlockRigidBody>().rigidBody);
        }
    }


    #region DetectCollision
    void DetectCollisions()
    {
        for(int i = 0;i < blockList.Count;i++)
            for(int j = i+1;j < blockList.Count;j++)
            {
				AddCollision(i, j);
            }
        
    }

    bool AddCollision(int index1,int index2)
    {
        //Fetching GameObjects and their rigidbodies
        GameObject r1 = blockList[index1];
        GameObject r2 = blockList[index2];
        RigidBodyClass rb1 = rigidBodies[index1];
        RigidBodyClass rb2 = rigidBodies[index2];

        //Getting normals and vertices of each object
        Vector3[] verti_a = new Vector3[8];
        Vector3[] verti_b = new Vector3[8];
        Vector3[] norma_a = new Vector3[6];
        Vector3[] norma_b = new Vector3[6];

        GetNormalAndVertices(r1, r2, ref norma_a, ref norma_b, ref verti_a, ref verti_b);

        float min_p = 10000;
        Vector3 min = new Vector3(0, 0, 0);
        bool facev = false;
        Vector3 ea = new Vector3(0, 0, 0);
        Vector3 eb = new Vector3(0, 0, 0);

        //normals of a -> 3
        int c = 0;
        for (int i = 0; i < norma_a.Length; i++)
        {
            if (norma_a[i][0] >= 0 && norma_a[i][1] >= 0 && norma_a[i][2] >= 0)
            {
                //Debug.Log(norma_a[i]);
                c++;
                float[] ans = IntersectsWhenProjected(verti_a, verti_b, norma_a[i]);

                if (ans[0] != 1)
                {
                    return false;
                }
                else
                {
                    if (ans[2] - ans[1] < min_p)
                    {
                        min_p = ans[2] - ans[1];
                        min = norma_a[i];
                        facev = true;
                        //print("FaceVertex");
                    }
                }
            }
        }

        //normals of b -> 3
        for (int j = 0; j < norma_b.Length; j++)
        {
            if (norma_b[j][0] >= 0 && norma_b[j][1] >= 0 && norma_b[j][2] >= 0)
            {
                //Debug.Log(norma_b[j]);
                c++;
                float[] ans = IntersectsWhenProjected(verti_a, verti_b, norma_b[j]);

                if (ans[0] != 1)
                {
                    return false;
                }
                else
                {
                    if (ans[2] - ans[1] < min_p)
                    {
                        min_p = ans[2] - ans[1];
                        min = norma_b[j];
                        facev = true;
                        //print("FaceVertex");
                    }
                    
                }
            }
        }

        //cross products of normals -> 9
        for (int i = 0; i < norma_a.Length; i++)
        {

            for (int j = 0; j < norma_b.Length; j++)
            {
                if (norma_a[i][0] >= 0 && norma_b[i][1] >= 0 && norma_b[i][2] >= 0)
                {
                    if (norma_b[j][0] >= 0 && norma_b[j][1] >= 0 && norma_b[j][2] >= 0)
                    {
                        c++;
                        float[] ans = IntersectsWhenProjected(verti_a, verti_b, Vector3.Cross(norma_a[i], norma_b[j]));

                        if (ans[0] != 1)
                        {
                            return false;
                        }
                        else
                        {
                            if (ans[2] - ans[1] < min_p)
                            {
                                min_p = ans[2] - ans[1];
                                min = norma_b[j];
                                facev = false;
                                ea = norma_a[i];
                                eb = norma_b[j];
                                //print("Edge/Edge");
                            }
                        }

                    }

                }

            }
        }

		Debug.Log("collision");
        //Vector3 temp = new Vector3((rb1.GetPosition().x - rb2.GetPosition().x) * min.x, (rb1.GetPosition().y - rb2.GetPosition().y) * min.y, (rb1.GetPosition().z - rb2.GetPosition().z) * min.z);
		Vector3 p1 = (rb2.GetPosition() + ((rb1.GetPosition() - rb2.GetPosition()) / 2));

		Vector3 temp = rb1.GetPosition() - rb2.GetPosition();

		//if (Vector3.Dot(temp, min) < 0)
		//{
		//	min = -1 * min;
		//}


        Collision col = new Collision(rb1, rb2, p1, min, ea, eb, facev);
        collisions.Add(col);
		
		

        //Collision Detected between i and j
        //rb1.CollisionTest();
        //rb2.CollisionTest();
        
        return true;
    }

    //Helper Functions below
    void GetNormalAndVertices(GameObject obj1,GameObject obj2,ref Vector3[] norma_a, ref Vector3[] norma_b, ref Vector3[] verti_a, ref Vector3[] verti_b)
    {
        //Getting Mesh from GameObjects
        Mesh a = obj1.GetComponent<MeshFilter>().mesh;
        Mesh b = obj2.GetComponent<MeshFilter>().mesh;

        //We just want the normals to the 6 edges
        norma_a[0] = a.normals[0];
        norma_a[1] = a.normals[4];
        norma_a[2] = a.normals[6];
        norma_a[3] = a.normals[12];
        norma_a[4] = a.normals[16];
        norma_a[5] = a.normals[20];

        norma_b[0] = b.normals[0];
        norma_b[1] = b.normals[4];
        norma_b[2] = b.normals[6];
        norma_b[3] = b.normals[12];
        norma_b[4] = b.normals[16];
        norma_b[5] = b.normals[20];

        for (int i = 0; i < 8; i++)
        {
            verti_a[i] = obj1.transform.TransformPoint(a.vertices[i]);
            verti_b[i] = obj2.transform.TransformPoint(b.vertices[i]);
        }
    }
    // aCorn and bCorn are arrays containing all corners (vertices) of the two OBBs
    private static float[] IntersectsWhenProjected(Vector3[] aCorn, Vector3[] bCorn, Vector3 axis)
    {


        float[] ans = new float[3];
        ans[0] = 1;
        ans[1] = 0;
        ans[2] = 100000;

        // Handles the cross product = {0,0,0} case
        if (axis == Vector3.zero)
            return ans;

        float aMin = float.MaxValue;
        float aMax = float.MinValue;
        float bMin = float.MaxValue;
        float bMax = float.MinValue;

        // Define two intervals, a and b. Calculate their min and max values
        for (int i = 0; i < 8; i++)
        {
            float aDist = Vector3.Dot(aCorn[i], axis);
            aMin = (aDist < aMin) ? aDist : aMin;
            aMax = (aDist > aMax) ? aDist : aMax;
            float bDist = Vector3.Dot(bCorn[i], axis);
            bMin = (bDist < bMin) ? bDist : bMin;
            bMax = (bDist > bMax) ? bDist : bMax;
        }

        // One-dimensional intersection test between a and b
        float longSpan = Mathf.Max(aMax, bMax) - Mathf.Min(aMin, bMin);
        float sumSpan = aMax - aMin + bMax - bMin;
        //return longSpan < sumSpan; // Change this to <= if you want the case were they are touching but not overlapping, to count as an intersection
        ans[0] = (longSpan <= sumSpan) ? (1) : (0);
        ans[1] = longSpan;
        ans[2] = sumSpan;
        return ans;
    }

    #endregion

    #region Colliding Contact
    void DefineCollisionType()
    {
        
        //Calculate which type of collision it is
        for (int i = 0;i < collisions.Count;i++)
            CollisionType(collisions[i]);
    }

    void CollidingContact()
    {
        double epsilon = 0.3;
        //For all colliding contacts...
        for (int i = 0; i < collidingContact.Count; i++)
        {
            Collision(collidingContact[i], epsilon);
        }
    }

    void Collision(Collision c,double epsilon)
    {
		Debug.Log("Performing colliding contact");
        Vector3 padot = VelocityAtPoint(c.obj1, c.CollisionLocation);
        Vector3 pbdot = VelocityAtPoint(c.obj2, c.CollisionLocation);
        Vector3 n = c.NormalDirection;
        Vector3 ra = c.CollisionLocation - c.obj1.GetPosition();
        Vector3 rb = c.CollisionLocation - c.obj2.GetPosition();

        double vrel = Vector3.Dot(n, (padot - pbdot));
        //Debug.Log(vrel);
        double numerator = -(1 + epsilon) * vrel;

        double term1 = c.obj1.GetMassInverse();
        double term2 = c.obj2.GetMassInverse(); 
        double term3 = Vector3.Dot(n,Vector3.Cross(c.obj1.GetIInverse().MultiplyVector((Vector3.Cross(ra, n))),ra));
        double term4 = Vector3.Dot(n, Vector3.Cross(c.obj2.GetIInverse().MultiplyVector((Vector3.Cross(rb, n))), rb));

        double j = numerator / (term1 + term2 + term3 + term4);
        Vector3 force = n * (float)j;
		//print(force);
		Debug.Log(n);
		Debug.Log(j);
		Debug.Log(force);


		c.obj1.AddLinearMomentum(force);
		c.obj2.AddLinearMomentum(-force);
	
		//c.obj1.AddForce(force);
		//c.obj2.AddForce(-force);

		//c.obj1.AddTorque(Vector3.Cross(ra, force));
		//c.obj2.AddTorque(-Vector3.Cross(ra, force));


		
        c.obj1.AddAngularMomentum(Vector3.Cross(ra, force));
        c.obj2.AddAngularMomentum(-Vector3.Cross(rb, force));
		/*
		c.obj1.RecomputeVelocity();
        c.obj2.RecomputeVelocity();
        */
	}

    //Helper Functions
    Vector3 VelocityAtPoint(RigidBodyClass rb,Vector3 point)
    {
        return (rb.GetVelocity() + Vector3.Cross(rb.GetAngularVelocity(), point - rb.GetPosition()));
    }
    void CollisionType(Collision c)
    {
        Vector3 padot = VelocityAtPoint(c.obj1, c.CollisionLocation);
        Vector3 pbdot = VelocityAtPoint(c.obj2, c.CollisionLocation);
        double vrel = Vector3.Dot(c.NormalDirection, (padot - pbdot));

        if (vrel > 5)
        {
            
        }
        //else if (vrel < -10)
        //{
			//restingContact.Add(c);
        //}
        else
        {
			Debug.Log("Colliding Contact");
            collidingContact.Add(c);
        }

    }

    #endregion

    #region RestingContact
    void RestingContact()
    {
        NMatrix AMatrix = new NMatrix(restingContact.Count, restingContact.Count);
        NVector fVector = new NVector(restingContact.Count);
        NVector bVector = new NVector(restingContact.Count);

        ComputeAMatrix(ref AMatrix);
        ComputeBVector(ref bVector);
        
        //Now I need a QP Solver. Not completely sure which one to use, and how. Uncomment when you do this
        //QP_Solve(AMatrix,bVector,ref fVector);

        for(int i = 0;i < restingContact.Count;i++)
        {
            double force = fVector.GetPosition(i);
            Vector3 normalDirection = restingContact[i].CollisionLocation;
            RigidBodyClass classA = restingContact[i].obj1;
            RigidBodyClass classB = restingContact[i].obj2;

            //Applying force in the positive direction
            classA.AddForce((float)force * normalDirection);
            Vector3 TorqueDistance = restingContact[i].CollisionLocation - classA.GetPosition();
            classA.AddTorque(Vector3.Scale(TorqueDistance,(float)force*normalDirection));

            force = -1 * force;
            //Applying force in the negative direction
            classB.AddForce((float)force * normalDirection);
            TorqueDistance = restingContact[i].CollisionLocation - classB.GetPosition();
            classB.AddTorque(Vector3.Scale(TorqueDistance, (float)force * normalDirection));

        }


    }

    //Helper Functions for Resting Contact
    void ComputeAMatrix(ref NMatrix aMatrix)
    {
        for(int i = 0;i < restingContact.Count;i++)
        {
            for (int j = 0; j < restingContact.Count; j++)
                aMatrix.SetPosition(i,j,ComputeAij(i,j));            
        }
    }
    float ComputeAij(int i,int j)
    {
        //Check if the bodies influence ea or not
        Collision ci = restingContact[i];
        Collision cj = restingContact[j];
        if ((ci.obj1 != cj.obj1) && (ci.obj1 != cj.obj2) && (ci.obj2 != cj.obj1) && (ci.obj2 != cj.obj2))
            return 0;
        else
        {
            RigidBodyClass a = ci.obj1;
            RigidBodyClass b = ci.obj2;
            Vector3 ni = ci.NormalDirection;
            Vector3 nj = cj.NormalDirection;
            Vector3 pi = ci.CollisionLocation;
            Vector3 pj = cj.CollisionLocation;
            Vector3 ra = pi - a.GetPosition();
            Vector3 rb = pi - a.GetPosition();

            Vector3 forceOnA = Vector3.zero;
            Vector3 torqueOnA = Vector3.zero;
            Vector3 forceOnB = Vector3.zero;
            Vector3 torqueOnB = Vector3.zero;

            if(ci.obj1 == cj.obj1)
            {
                forceOnA = nj;
                torqueOnA = Vector3.Cross((pj - a.GetPosition()), nj);
            }
            else if(ci.obj1 == cj.obj2)
            {
                forceOnA = -nj;
                torqueOnA = Vector3.Cross((pj - a.GetPosition()), nj);
            }

            if (ci.obj2 == cj.obj1)
            {
                forceOnB = nj;
                torqueOnB = Vector3.Cross((pj - b.GetPosition()), nj);
            }
            else if (ci.obj2 == cj.obj2)
            {
                forceOnB = -nj;
                torqueOnB = Vector3.Cross((pj - b.GetPosition()), nj);
            }

            Vector3 a_linear = forceOnA / (float)a.GetMass();
            Vector3 a_angular = Vector3.Cross(a.GetIInverse().MultiplyVector(torqueOnA),ra);
            Vector3 b_linear = forceOnB / (float)b.GetMass();
            Vector3 b_angular = Vector3.Cross(b.GetIInverse().MultiplyVector(torqueOnB), rb);

            float result = Vector3.Dot(ni, ((a_linear + a_angular) - (b_linear + b_angular)));

            return result;
        }
    }


    void ComputeBVector(ref NVector bVector)
    {

        for(int i = 0;i < restingContact.Count;i++)
        {
            Collision c = restingContact[i];
            RigidBodyClass obj1 = c.obj1;
            RigidBodyClass obj2 = c.obj2;
            Vector3 normal = c.NormalDirection;
            Vector3 ADistanceFromCenter = c.CollisionLocation - obj1.GetPosition();
            Vector3 BDistanceFromCenter = c.CollisionLocation - obj2.GetPosition();

            Vector3 AExtForce = obj1.GetForce();
            Vector3 BExtForce = obj2.GetForce();
            Vector3 AExtTorque = obj1.GetTorque();
            Vector3 BExtTorque = obj2.GetTorque();

            Vector3 AExtPart = AExtForce / (float)obj1.GetMass() + (Vector3.Cross((obj1.GetIInverse().MultiplyVector(AExtTorque)),ADistanceFromCenter));
            Vector3 BExtPart = BExtForce / (float)obj2.GetMass() + (Vector3.Cross((obj2.GetIInverse().MultiplyVector(BExtTorque)), BDistanceFromCenter));

            Vector3 AVelPart = (Vector3.Cross(obj1.GetAngularVelocity(), Vector3.Cross(obj1.GetAngularVelocity(), ADistanceFromCenter)));
            AVelPart += Vector3.Cross(obj1.GetIInverse().MultiplyVector(Vector3.Cross(obj1.GetAngularMomentum(),obj1.GetAngularVelocity())), ADistanceFromCenter);

            Vector3 BVelPart = (Vector3.Cross(obj2.GetAngularVelocity(), Vector3.Cross(obj2.GetAngularVelocity(), BDistanceFromCenter)));
            BVelPart += Vector3.Cross(obj2.GetIInverse().MultiplyVector(Vector3.Cross(obj2.GetAngularMomentum(), obj2.GetAngularVelocity())), BDistanceFromCenter);

            Vector3 PositionPartCombination = ((AExtPart + AVelPart) - (BExtPart - BVelPart));
            double k1 = Vector3.Dot(normal, PositionPartCombination);
            Vector3 ndot = ComputeNDot(c);
            double k2 = 2 * Vector3.Dot(ndot, (VelocityOfPointOnBody(obj1, c.CollisionLocation))-(VelocityOfPointOnBody(obj2,c.CollisionLocation)));
            bVector.SetPosition(i, (k1 + k2));
        }

        

    }
    Vector3 VelocityOfPointOnBody(RigidBodyClass rb, Vector3 point)
    {
        return rb.GetVelocity() + Vector3.Cross(rb.GetAngularVelocity(),(point - rb.GetPosition()));
    }
    Vector3 ComputeNDot(Collision collision)
    {
        Vector3 result = new Vector3(0, 0, 0);
        if(collision.IsVertexFaceContact)
        {
            result = Vector3.Cross(collision.obj2.GetAngularVelocity(), collision.NormalDirection);

        }
        else
        {
            Vector3 edgeAdot = Vector3.Cross(collision.obj1.GetAngularVelocity(), collision.Edge1Direction);
            Vector3 edgeBdot = Vector3.Cross(collision.obj2.GetAngularVelocity(), collision.Edge2Direction);
            Vector3 normal1 = Vector3.Cross(collision.Edge1Direction, collision.Edge2Direction);
            Vector3 ZDerivative = Vector3.Cross(edgeAdot, collision.Edge2Direction) + Vector3.Cross(collision.Edge1Direction, edgeBdot);
            double normalLength = normal1.magnitude;
            normal1.Normalize();
            result = (ZDerivative - Vector3.Cross(Vector3.Cross(ZDerivative, normal1),normal1));
            result = result / (float)normalLength;
            
        }
        return result;
    }

    void QP_Solve(NMatrix AMatrix,NVector bVector, ref NVector fVector)
    {
        //I'm unable to find a good way of solving QP problems for now. I'm going to assume that since 
    }
    #endregion
}
