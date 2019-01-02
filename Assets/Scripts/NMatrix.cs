using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NMatrix {

    float[,] matrix;

    public NMatrix(int size1,int size2)
    {
        matrix = new float[size1, size2];
    }

    public void SetPosition(int loc1,int loc2,float value) { matrix[loc1, loc2] = value; }

    public float GetPosition(int loc1,int loc2) { return matrix[loc1, loc2]; }
}
