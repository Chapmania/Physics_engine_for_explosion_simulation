using System.Collections;
using System.Collections.Generic;
using UnityEngine;


public class BlockRigidBody : MonoBehaviour {

    public RigidBodyClass rigidBody = new RigidBodyClass();

    [SerializeField]
    bool IsFloor = false;

	[SerializeField]
	bool IsStone = false;

    //Initializing all Derived Values and constant Values
    void Awake()
    {
        //Setting rigidBody initial Values
        CalculateIbody();
        rigidBody.SetPosition(gameObject.transform.position);
        rigidBody.SetRotation(gameObject.transform.rotation);
		if (!IsFloor)
        { //rigidBody.AddForce(new Vector3(0, -9.81f, 0));
            rigidBody.SetMassInverse(1 / rigidBody.GetMass());
        }
        else
        {
            rigidBody.SetMassInverse(0);
            rigidBody.SetIbodyInverse(Matrix4x4.zero);
        }

		if (IsStone)
		{
			rigidBody.SetMass(1);
			rigidBody.AddForce(new Vector3(15, 0, 0));
		}

    }

    void Start()
    {
        
    }

    //Happens every frame
    void Update()
    {
        gameObject.transform.position = rigidBody.CalculateLinearMotion();
        gameObject.transform.rotation = rigidBody.CalculateAngularMotion();
    }
 
    //Generating IBody Values for the block
    void CalculateIbody()
    {

        Collider m_Collider = GetComponent<Collider>();
        Vector3 size = m_Collider.bounds.size;

        //Calculate the Ibody of a block.
        Matrix4x4 Ibody = Matrix4x4.identity;
        Ibody[0, 0] = (float)((size.y * size.y + size.z * size.z) * rigidBody.GetMass() / 12);
        Ibody[1, 1] = (float)((size.x * size.x + size.z * size.z) * rigidBody.GetMass() / 12);
        Ibody[2, 2] = (float)((size.y * size.y + size.x * size.x) * rigidBody.GetMass() / 12);

        rigidBody.setIBody(Ibody);

    }


}
