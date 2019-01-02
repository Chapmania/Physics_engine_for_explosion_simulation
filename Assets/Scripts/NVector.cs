using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NVector{

    float[] vector;

    public NVector(int size)
    {
        vector = new float[size];
    }

    public double GetPosition(int index)    { return vector[index]; }

    public void SetPosition(int index, double value) { vector[index] = (float)value; }

    public void SetZero()
    {
        for (int i = 0; i < vector.Length; i++)
            vector[i] = 0;
    }
}
