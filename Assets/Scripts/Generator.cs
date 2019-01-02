using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Generator : MonoBehaviour {


	// Use this for initialization
	void Start () {
        
	}
	
    void RandomForce()
    {
        GameObject[] listofBlocks = GameObject.FindGameObjectsWithTag("Block");
        for (int i = 0; i < 10; i++)
        {
            int randBlock = Random.Range(0, listofBlocks.Length);
            print(randBlock);
            RigidBodyClass rb = listofBlocks[i].GetComponent<BlockRigidBody>().rigidBody;
            Vector3 randomForce = new Vector3(Random.value, Random.value, 0);
            rb.AddForce(randomForce);
        }
    }
	// Update is called once per frame
	void Update () {
		
	}
}
