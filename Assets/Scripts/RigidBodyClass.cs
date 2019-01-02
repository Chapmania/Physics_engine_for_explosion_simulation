using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RigidBodyClass {

    //Constants
    double Mass = 5;
    double MassInverse;
    Matrix4x4 Ibody;
    Matrix4x4 IbodyInverse;

    //State Variables
    Vector3 Position;
    Quaternion Rotation;
    Vector3 LinearMomentum = new Vector3(0, 0, 0);
    Vector3 AngularMomentum = new Vector3(0, 0, 0);

    //Derived Variables
    Matrix4x4 InertiaInverse;
    Vector3 Velocity;
    Vector3 AngularVelocity;

    //Computed Variables
    Vector3 Torque = new Vector3(0, 0, 0);
    Vector3 Force = new Vector3(0, 0, 0);

    //Setters
    public void setIBody(Matrix4x4 IBody) { Ibody = IBody; IbodyInverse = Ibody.inverse; }
    public void SetPosition(Vector3 position) { Position = position; }
    public void SetRotation(Quaternion rotation) { Rotation = rotation; }
    public void SetMass(double mass) { Mass = mass; }
    public void SetMassInverse(double massInverse) { MassInverse = massInverse; }
    public void SetIbodyInverse(Matrix4x4 Ibodyvalue) { IbodyInverse = Ibodyvalue; }
    public void AddAngularMomentum(Vector3 input) { AngularMomentum += input; }
    public void AddLinearMomentum(Vector3 input ) { LinearMomentum += input; }
    public void AddForce(Vector3 newForce) { Force += newForce; }
    public void AddTorque(Vector3 newTorque) { Torque += newTorque; }

    //Getters
    public double GetMass() { return Mass; }
    public Vector3 GetPosition() { return Position; }
    public Vector3 GetAngularVelocity() { return AngularVelocity; }
    public Vector3 GetForce() { return Force; }
    public Vector3 GetTorque() { return Torque; }
    public Matrix4x4 GetIInverse() { return InertiaInverse; }
    public Vector3 GetAngularMomentum() { return AngularMomentum; }
    public Vector3 GetVelocity() { return Velocity;  }
    public double GetMassInverse() { return MassInverse; }


    //Linear and Angular Motion Calculations
    public Vector3 CalculateLinearMotion()
    {
        //Update Linear Momentum using Force
        LinearMomentum = LinearMomentum + Time.deltaTime * Force;

        //Calculate Velocity using Linear Momentum
        Velocity = LinearMomentum * (float)MassInverse;

        //Calculate New Position using Velocity
        Position = Position + Velocity * Time.deltaTime;
        return Position;
    }
    public Quaternion CalculateAngularMotion()
    {
        //Update Angular Momentum using Torque
        AngularMomentum = AngularMomentum + Time.deltaTime * Torque;

        //Calculating InverseInertia using Angular Momentum
        Matrix4x4 RotationMatrix = Matrix4x4.Rotate(Rotation);
        InertiaInverse = RotationMatrix * IbodyInverse * Matrix4x4.Transpose(RotationMatrix);
        
        //Calculate Angular Velocity using InverseInertia
        AngularVelocity = InertiaInverse * AngularMomentum;
        //print(AngularVelocity);

        //Calculate Change in Rotation using Angular Velocity
        Matrix4x4 ChangeInRotation = Vector3ToMatrix(AngularVelocity) * RotationMatrix;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                RotationMatrix[i, j] = RotationMatrix[i, j] + (ChangeInRotation[i, j] * Time.deltaTime);
            }
        RotationMatrix[3, 3] = 1;
        Rotation = RotationMatrix.rotation;
        return Rotation;
    }

    //Helper Function : Convert a Vector3 to Matrix4x4
    Matrix4x4 Vector3ToMatrix(Vector3 vector3)
    {
        Matrix4x4 result = Matrix4x4.identity;
        result[0, 0] = 0;
        result[0, 1] = -vector3[2];
        result[0, 2] = vector3[1];
        result[1, 0] = vector3[2];
        result[1, 1] = 0;
        result[1, 2] = -vector3[0];
        result[2, 0] = -vector3[1];
        result[2, 1] = vector3[0];
        result[2, 2] = 0;
        return result;
    }

    public void RecomputeVelocity()
    {
        Velocity = LinearMomentum / (float)Mass;
        //Calculating InverseInertia using Angular Momentum
        Matrix4x4 RotationMatrix = Matrix4x4.Rotate(Rotation);
        InertiaInverse = RotationMatrix * IbodyInverse * Matrix4x4.Transpose(RotationMatrix);

        //Calculate Angular Velocity using InverseInertia
        AngularVelocity = InertiaInverse * AngularMomentum;
    }
    public void CollisionTest()
    {
        LinearMomentum = Vector3.zero;
        Force = Vector3.zero;
        AngularMomentum = Vector3.zero;
        Torque = Vector3.zero;
    }

}
