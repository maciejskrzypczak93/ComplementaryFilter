package com.example.marek.nazwa;

import android.hardware.SensorManager;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.support.v7.widget.Toolbar;
import android.util.Log;
import android.view.View;
import android.view.Menu;
import android.view.MenuItem;
import android.widget.Button;
import android.widget.TextView;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;

public class MainActivity extends AppCompatActivity {

    ArrayList<Float> t_acc = new ArrayList<Float>();
    ArrayList<Integer> t_acc_int = new ArrayList<Integer>();
    ArrayList<Float> x_acc = new ArrayList<Float>();
    ArrayList<Float> y_acc = new ArrayList<Float>();
    ArrayList<Float> z_acc = new ArrayList<Float>();
    ArrayList<Float> t_gyro = new ArrayList<Float>();
    ArrayList<Integer> t_gyro_int = new ArrayList<Integer>();
    ArrayList<Float> x_gyro = new ArrayList<Float>();
    ArrayList<Float> y_gyro = new ArrayList<Float>();
    ArrayList<Float> z_gyro = new ArrayList<Float>();
    ArrayList<Float> Pitch_orient = new ArrayList<Float>();
    ArrayList<Float> Roll_orient = new ArrayList<Float>();
    ArrayList<Float> Yaw_orient = new ArrayList<Float>();
    ArrayList<Float> Pitch_gyro = new ArrayList<Float>();
    ArrayList<Float> Roll_gyro = new ArrayList<Float>();
    ArrayList<Float> Yaw_gyro = new ArrayList<Float>();
    ArrayList<Float> Pitch = new ArrayList<Float>();
    ArrayList<Float> Roll = new ArrayList<Float>();
    ArrayList<Float> Yaw = new ArrayList<Float>();

    ArrayList<Float> t_orient = new ArrayList<Float>();
    ArrayList<Integer> t_orient_int = new ArrayList<Integer>();
    ArrayList<Float> x_orient = new ArrayList<Float>();
    ArrayList<Float> y_orient = new ArrayList<Float>();
    ArrayList<Float> z_orient = new ArrayList<Float>();

    ArrayList<float[]> quaternion_acm = new ArrayList<>();
    ArrayList<float[]> quaternion_acm1 = new ArrayList<>();
    ArrayList<Float> t_data = new ArrayList<Float>();
    ArrayList<Float> w_data = new ArrayList<Float>();
    ArrayList<Float> x_data = new ArrayList<Float>();
    ArrayList<Float> y_data = new ArrayList<Float>();
    ArrayList<Float> z_data = new ArrayList<Float>();
    ArrayList<Integer> t_data_int = new ArrayList<Integer>();
    ArrayList<Float> V_Length = new ArrayList<Float>();
    ArrayList<Float> V_Length_sort = new ArrayList<Float>();
    float[] a =null;
    float[] b =null;
    float[] c =null;
    float[] d =null;
    float[] w =null;
    float[] x =null;
    float[] y =null;
    float[] z =null;
    float[] QuatGyroW =null;
    float[] QuatGyroX =null;
    float[] QuatGyroY =null;
    float[] QuatGyroZ =null;
    ArrayList<float[]> quaternion_gyro = new ArrayList<>();
    ArrayList<float[]> quaternion_gyro_mul = new ArrayList<>();
    ArrayList<float[]> quaternion_gyro1 = new ArrayList<>();
    ArrayList<float[]> quaternion_gyro_mul1 = new ArrayList<>();
    ArrayList<float[]> quaternion_gyro_logmap = new ArrayList<>();
    ArrayList<float[]> quaternion_gyro_logmap1 = new ArrayList<>();
    ArrayList<float[]> quaternion_orient_logmap = new ArrayList<>();
    ArrayList<float[]> quaternion_orient_logmap1 = new ArrayList<>();
    ArrayList<float[]> vector_filter = new ArrayList<>();
    ArrayList<float[]> vector_filter1 = new ArrayList<>();
    ArrayList<float[]> quaternion_exp = new ArrayList<>();
    ArrayList<float[]> quaternion_exp1 = new ArrayList<>();
    ArrayList<float[]> QxR = new ArrayList<>();
    ArrayList<float[]> QxR1 = new ArrayList<>();
    ArrayList<float[]> InvSense = new ArrayList<>();
    ArrayList<float[]> InvSense1 = new ArrayList<>();
    ArrayList<float[]> QxRxInvSense = new ArrayList<>();
    ArrayList<float[]> QxRxInvSense1 = new ArrayList<>();
    ArrayList<float[]> QxRxInvSenseLog = new ArrayList<>();
    ArrayList<float[]> QxRxInvSenseLog1 = new ArrayList<>();
    ArrayList<float[]> Sense = new ArrayList<>();
    ArrayList<float[]> Sense1 = new ArrayList<>();

    float[] OrientQuat =null;
    float[] qacm = new float [3];
    float[] qgyro = new float [3];
    float[] GyroQuat =new float [4];
    float[] w_filter =null;
    float[] x_filter =null;
    float[] y_filter =null;
    float[] z_filter =null;
    float[] w_filter_exp =null;
    float[] x_filter_exp =null;
    float[] y_filter_exp =null;
    float[] z_filter_exp =null;
    float[] w_QxR =null;
    float[] x_QxR =null;
    float[] y_QxR =null;
    float[] z_QxR =null;
    float[] w_QxData =null;
    float[] x_QxData =null;
    float[] y_QxData =null;
    float[] z_QxData =null;
    float[] w_LogQxData =null;
    float[] x_LogQxData =null;
    float[] y_LogQxData =null;
    float[] z_LogQxData =null;
    int k=0;
    float error=0;

    int licznik_przycisk=0;
    Float gyro_sum_x= Float.valueOf(0);
    Float gyro_sum_y=Float.valueOf(0);
    Float gyro_sum_z=Float.valueOf(0);
    static final double ms2s = 1e-3;
    double PI = 3.14159265359;
    private static final float eps = 1e-9f;
    float waga=0.98f;

    float dt;

    //odczyt danych
    void filereader(java.lang.String filename, String device){
        BufferedReader reader = null;
        String[] line_string;
        String mLine=null;
        float[] tab;
        try {
            reader = new BufferedReader(
                        new InputStreamReader(getAssets().open(filename)));
            while ((mLine = reader.readLine()) != null) {
                //process line
                line_string = mLine.split(" ");
                tab = new float[line_string.length];

                for (int i=0; i<line_string.length; i++){
                    tab[i]=Float.valueOf(line_string[i]);
                }

                if(device=="acc"){
                    t_acc.add(tab[0]);
                    t_acc_int.add((int)tab[0]);
                    x_acc.add(tab[1]);
                    y_acc.add(tab[2]);
                    z_acc.add(tab[3]);
                }
                if(device=="gyro"){
                    t_gyro.add(tab[0]);
                    t_gyro_int.add((int)tab[0]);
                    x_gyro.add(tab[1]);
                    y_gyro.add(tab[2]);
                    z_gyro.add(tab[3]);
                }
                if(device=="orient"){
                    t_orient.add(tab[0]);
                    t_orient_int.add((int)tab[0]);
                    x_orient.add(tab[1]);
                    y_orient.add(tab[2]);
                    z_orient.add(tab[3]);
                }
                if(device=="data"){
                    t_data.add(tab[0]);
                    t_data_int.add((int)tab[0]);
                    w_data.add(tab[4]);
                    x_data.add(tab[5]);
                    y_data.add(tab[6]);
                    z_data.add(tab[7]);
                }
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static float[] quatMul(float[] q1, float[] q2){
        float[] ret = new float[4];
        ret[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
        ret[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2];
        ret[2] = q1[0]*q2[2] + q1[2]*q2[0] + q1[3]*q2[1] - q1[1]*q2[3];
        ret[3] = q1[0]*q2[3] + q1[3]*q2[0] + q1[1]*q2[2] - q1[2]*q2[1];
        return ret;
    }

    private static float[] quatInv(float[] q){
        float[] ret = new float[4];
        ret[0] = q[0];
        ret[1] = -q[1];
        ret[2] = -q[2];
        ret[3] = -q[3];

        return ret;
    }


    void compQuatFromGyro(float[] gyroVals,
                          float dt,
                          float[] quat)
    {
        float[] gyroValsNorm = new float[]{0.0f, 0.0f, 0.0f};
        float norm = (float)Math.sqrt(gyroVals[0] * gyroVals[0] +
                gyroVals[1] * gyroVals[1] +
                gyroVals[2] * gyroVals[2]);
        if(norm > eps){
            for(int v = 0; v < 3; v++){
                gyroValsNorm[v] = gyroVals[v] / norm;
            }
        }

        float thetaOverTwo = norm * dt / 2.0f;
        float sinThetaOverTwo = (float)Math.sin(thetaOverTwo);
        float cosThetaOverTwo = (float)Math.cos(thetaOverTwo);

        quat[0] = cosThetaOverTwo;
        quat[1] = sinThetaOverTwo * gyroValsNorm[0];
        quat[2] = sinThetaOverTwo * gyroValsNorm[1];
        quat[3] = sinThetaOverTwo * gyroValsNorm[2];
    }

    private static float[] logMap(float[] quat){
        float[] res = new float[3];
        double qvNorm = Math.sqrt
                ((quat[1]*quat[1]+quat[2]*quat[2]+quat[3]*quat[3]));
        if(qvNorm > 1e-6){
            res[0] = quat[0]/(float)qvNorm;
            res[1] = quat[1]/(float)qvNorm;
            res[2] = quat[2]/(float)qvNorm;
        }
        else{
            // 1/sqrt(3), so norm = 1
            res[0] = 0.57735026919f;
            res[1] = 0.57735026919f;
            res[2] = 0.57735026919f;
        }
        double acosQw = Math.acos(quat[0]);
        res[0] *= 2.0*acosQw;
        res[1] *= 2.0*acosQw;
        res[2] *= 2.0*acosQw;
        return res;
    }


    private static float[] expMap(float[] vec){
        double arg = 0.5*Math.sqrt
                (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
        double sincArg = 1.0;
        if(arg>1e-6){
            sincArg = Math.sin(arg)/arg;
        }
        else{
            sincArg = 1 - arg*arg/6 + Math.pow(arg, 4)/120;
        }
        double cosArg = Math.cos(arg);
        float[] ret = new float[4];
        ret[0] = (float)cosArg;
        ret[1] = 0.5f*(float)sincArg*vec[0];
        ret[2] = 0.5f*(float)sincArg*vec[1];
        ret[3] = 0.5f*(float)sincArg*vec[2];

        //normalize and unify
        float n = 1.0f / (float)Math.sqrt
                (ret[0]*ret[0]+ret[1]*ret[1]+ret[2]*ret[2]+ret[3]*ret[3]);
        ret[0] *= n;
        ret[1] *= n;
        ret[2] *= n;
        ret[3] *= n;
        double eps = 1e-9;
        if(ret[0]<0.0 ||
                Math.abs(ret[0])<eps && ret[3] <0.0 ||
                Math.abs(ret[0])<eps && Math.abs(ret[3])<eps && ret[2] <0.0 ||
                Math.abs(ret[0])<eps && Math.abs(ret[3])<eps &&
                        Math.abs(ret[2])<eps && ret[1] <0.0)
        {
            ret[0] = -ret[0];
            ret[1] = -ret[1];
            ret[2] = -ret[2];
            ret[3] = -ret[3];
        }
        //////
        return ret;
    }


    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
        Toolbar toolbar = (Toolbar) findViewById(R.id.toolbar);
        setSupportActionBar(toolbar);


        BufferedReader reader = null;
        BufferedReader reader1 = null;
        BufferedReader reader2 = null;
        BufferedReader reader3 = null;
        BufferedReader reader4 = null;
        String mLine=null;
        String mLine1=null;
        String mLine2=null;
        String mLine3=null;
        String mLine4=null;


        filereader("acc.log","acc");
        filereader("gyro.log","gyro");
        filereader("orient.log","orient");
        filereader("Xsense.data","data");

        for (int i=0; i<t_gyro_int.size();i++){  //t_gyro-int.size();
            //tworzenie wektora xyz
            qgyro[0]=(float)x_gyro.get(i);
            qgyro[1]=(float)y_gyro.get(i);
            qgyro[2]=(float)z_gyro.get(i);

            //dt dla pierwszej probki
            if (i==0) {
                dt=1*(float)ms2s;
            }
            //dt dla pozostalych probek
            else{
                dt = (t_gyro_int.get(i) - t_gyro_int.get(i-1))*(float)ms2s;
            }
            //obliczenie kwaternionu
            compQuatFromGyro(qgyro, dt , GyroQuat);

            //zapis do arraylisty
            if(quaternion_gyro1.isEmpty()){
                quaternion_gyro1.add(GyroQuat);
                quaternion_gyro.add(quaternion_gyro1.get(i).clone());
            }
            else {
                quaternion_gyro1.add(quaternion_gyro1.size() - 1, GyroQuat);
                quaternion_gyro.add(quaternion_gyro1.get(i).clone());
            }
        }

        //mnozenie kwaternionow poprzedni razy nastepny

        QuatGyroW = new float[quaternion_gyro.size()];
        QuatGyroX = new float[quaternion_gyro.size()];
        QuatGyroY = new float[quaternion_gyro.size()];
        QuatGyroZ = new float[quaternion_gyro.size()];
        for(int i=0;i<quaternion_gyro.size();i++) {

            if (i == 0) {
                quaternion_gyro_mul1.add(quaternion_gyro.get(i));
                quaternion_gyro_mul.add(quaternion_gyro_mul1.get(i).clone());
            }
            if (i > 0) {
                quaternion_gyro_mul1.add(quatMul(quaternion_gyro.get(i - 1), quaternion_gyro.get(i)));
                quaternion_gyro_mul.add(quaternion_gyro_mul1.get(i).clone());
                //////////////////////////////////////////////// potrzebne do wyswietlania
                QuatGyroW[i] = quaternion_gyro_mul.get(i)[0];
                QuatGyroX[i] = quaternion_gyro_mul.get(i)[1];
                QuatGyroY[i] = quaternion_gyro_mul.get(i)[2];
                QuatGyroZ[i] = quaternion_gyro_mul.get(i)[3];
                ////////////////////////////////////////////////
            }
        }
        //LogMap gyro
        /////////////////////////
        w = new float[quaternion_gyro.size()];
        x = new float[quaternion_gyro.size()];
        y = new float[quaternion_gyro.size()];
        z = new float[quaternion_gyro.size()];
        ////////////////////////
        for(int i=0;i<quaternion_gyro.size();i++) {
            quaternion_gyro_logmap1.add(logMap(quaternion_gyro_mul.get(i)));
            quaternion_gyro_logmap.add(quaternion_gyro_logmap1.get(i).clone());
            ////////////////////////////////////////////
            w[i] = 0;
            x[i] = quaternion_gyro_logmap.get(i)[0];
            y[i] = quaternion_gyro_logmap.get(i)[1];
            z[i] = quaternion_gyro_logmap.get(i)[2];
            ////////////////////////////////////////////
        }

        ///////////////////////////////////////////////////////////

        for(int i=0;i<t_orient.size();i++){  //t_orient.size()
            //tworzenie wektora xyz
            qacm[0]=(float)x_orient.get(i);
            qacm[1]=(float)y_orient.get(i);
            qacm[2]=(float)z_orient.get(i);
            //inicjalizacja wektora
            if(OrientQuat==null)
            {
                OrientQuat = new float[4];
            }
            //obliczenie kwaternionu
            SensorManager.getQuaternionFromVector(OrientQuat,qacm);
            //zapis do arraylisty
            quaternion_acm1.add(OrientQuat);
            quaternion_acm.add(quaternion_acm1.get
                    (quaternion_acm1.size()-1).clone());

        }
        //////////////////////////////
        a = new float[quaternion_acm.size()];
        b = new float[quaternion_acm.size()];
        c = new float[quaternion_acm.size()];
        d = new float[quaternion_acm.size()];
        //////////////////////////////
        //LogMap orient
        for(int i=0;i<quaternion_acm.size();i++){
            quaternion_orient_logmap1.add(logMap(quaternion_acm.get(i)));
            quaternion_orient_logmap.add(quaternion_orient_logmap1.get(i).clone());
            ////////////////////////////////
            a[i]=0;
            b[i]=quaternion_gyro_logmap.get(i)[0];
            c[i]=quaternion_gyro_logmap.get(i)[1];
            d[i]=quaternion_gyro_logmap.get(i)[2];
            /////////////////////////////////
        }



        float[] xyz_filter = new float[3];
        int tSize=0;
        int i_orient=0;
        int i_gyro=0;
        if(t_orient_int.get(t_orient_int.size()-1)>=t_gyro_int.get(t_gyro_int.size()-1)){
            tSize=(t_orient_int.get(t_orient_int.size()-1));
        }
        else
            tSize=t_gyro_int.get(t_gyro_int.size()-1);

        w_filter = new float[tSize];
        x_filter = new float[tSize];
        y_filter = new float[tSize];
        z_filter = new float[tSize];
        w_filter_exp = new float[tSize];
        x_filter_exp = new float[tSize];
        y_filter_exp = new float[tSize];
        z_filter_exp = new float[tSize];
        for (int i=0; i<tSize-1 ;i++) {
            if(i_gyro>=quaternion_gyro_logmap.size()-1){
                i_gyro=quaternion_gyro_logmap.size()-1;
            }
            if(i_gyro>=quaternion_orient_logmap.size()-1){
                i_gyro=quaternion_orient_logmap.size()-1;
            }
            ///mnozenie z waga
            w_filter[i] = 0;
            x_filter[i] = waga * quaternion_gyro_logmap.get(i_gyro)[0] +
                    (1 - waga) * quaternion_orient_logmap.get(i_orient)[0];
            y_filter[i] = waga * quaternion_gyro_logmap.get(i_gyro)[1] +
                    (1 - waga) * quaternion_orient_logmap.get(i_orient)[1];
            z_filter[i] = waga * quaternion_gyro_logmap.get(i_gyro)[2] +
                    (1 - waga) * quaternion_orient_logmap.get(i_orient)[2];
            ////////////////////////////////
            xyz_filter[0]=waga * quaternion_gyro_logmap.get(i_gyro)[0] +
                    (1 - waga) * quaternion_orient_logmap.get(i_orient)[0];
            xyz_filter[1]=waga * quaternion_gyro_logmap.get(i_gyro)[1] +
                    (1 - waga) * quaternion_orient_logmap.get(i_orient)[1];
            xyz_filter[2]=waga * quaternion_gyro_logmap.get(i_gyro)[2] +
                    (1 - waga) * quaternion_orient_logmap.get(i_orient)[2];
            //////////////////////////////
            if(i<t_orient_int.size()-1){
                if(t_orient_int.get(i+1)==t_orient_int.get(i)){
                    i_orient++;
                }
            }
            if(i<t_gyro_int.size()-1){
                if(t_gyro_int.get(i+1)==t_gyro_int.get(i)){
                    i_gyro++;
                }
            }
            if (t_orient_int.contains(i)&&t_gyro_int.contains(i)){
                i_orient++;
                i_gyro++;
            }
            else if (t_orient_int.contains(i)&&t_gyro_int.contains(i)==false){
                i_orient++;
            }
            else if (t_orient_int.contains(i)==false&&t_gyro_int.contains(i)){
                i_gyro++;
            }

            vector_filter1.add(xyz_filter);
            vector_filter.add(vector_filter1.get(i).clone());
            //exp
            quaternion_exp1.add(expMap(vector_filter.get(i)));
            quaternion_exp.add(quaternion_exp1.get(i).clone());
            /////////////////
            w_filter_exp[i] = quaternion_exp.get(i)[0];
            x_filter_exp[i] = quaternion_exp.get(i)[1];
            y_filter_exp[i] = quaternion_exp.get(i)[2];
            z_filter_exp[i] = quaternion_exp.get(i)[3];
            //////////////////////////
        }
        //rotation matrix to quaternion
        int m00, m11, m22, m21, m12, m02, m20, m10, m01;
        float qw, qx, qy, qz;
        m00 = 0;
        m01 = 1;
        m02 = 0;
        m10 = 1;
        m11 = 0;
        m12 = 0;
        m20 = 0;
        m21 = 0;
        m22 = -1;
        float tr = m00 + m11 + m22;

        if (tr > 0) {
            float S = (float)Math.sqrt(tr+1.0) * 2; // S=4*qw
            qw = 0.25f * S;
            qx = (m21 - m12) / S;
            qy = (m02 - m20) / S;
            qz = (m10 - m01) / S;
        } else if ((m00 > m11)&(m00 > m22)) {
            float S = (float)Math.sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx
            qw = (m21 - m12) / S;
            qx = 0.25f * S;
            qy = (m01 + m10) / S;
            qz = (m02 + m20) / S;
        } else if (m11 > m22) {
            float S = (float)Math.sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy
            qw = (m02 - m20) / S;
            qx = (m01 + m10) / S;
            qy = 0.25f * S;
            qz = (m12 + m21) / S;
        } else {
            float S = (float)Math.sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz
            qw = (m10 - m01) / S;
            qx = (m02 + m20) / S;
            qy = (m12 + m21) / S;
            qz = 0.25f * S;
        }
        float[] matrixQuat = new float[4];
        matrixQuat[0] = qw;
        matrixQuat[1] = qx;
        matrixQuat[2] = qy;
        matrixQuat[3] = qz;

        /////////////////////////////

        //Q_filter*Rotation matrix Quaternion
        for(int i=0; i<quaternion_exp.size(); i++){
            QxR1.add(quatMul(matrixQuat, quaternion_exp.get(i)));
            QxR.add(QxR1.get(i).clone());
        }
        //////////////////////////////////
        w_QxR = new float[QxR.size()];
        x_QxR = new float[QxR.size()];
        y_QxR = new float[QxR.size()];
        z_QxR = new float[QxR.size()];
        /////////////////////////////////////////////
        for(int i=0; i<QxR.size();i++){

            w_QxR[i]=QxR.get(i)[0];
            x_QxR[i]=QxR.get(i)[1];
            y_QxR[i]=QxR.get(i)[2];
            z_QxR[i]=QxR.get(i)[3];
        }
        int size=0;
        int i_data=0;
        float[] wektor = new float[4];
        if(t_orient_int.get(t_orient_int.size()-1)>=t_gyro_int.get(t_gyro_int.size()-1)){
            size=t_orient_int.get(t_orient_int.size()-1);
        }
        else
            size=t_gyro_int.get(t_gyro_int.size()-1);
        for(int i=0; i<size-1; i++){
            if(i==t_data_int.get(i_data)){
                i_data++;
            }
            wektor[0]= w_data.get(i_data);
            wektor[1]= x_data.get(i_data);
            wektor[2]= y_data.get(i_data);
            wektor[3]= z_data.get(i_data);
            Sense1.add(wektor);
            Sense.add(Sense1.get(i).clone());
        }


        float[] invSense = new float[4];
        //invert quaternion data (Xsense)
        for (int i=0; i<Sense.size();i++) {
            invSense[0] = 0;
            invSense[1] = -Sense.get(i)[1];
            invSense[2] = -Sense.get(i)[2];
            invSense[3] = -Sense.get(i)[3];
            InvSense1.add(invSense);
            InvSense.add(InvSense1.get(i).clone());
        }
        //invert quaternion data (Xsense)


        w_QxData = new float[InvSense.size()];
        x_QxData = new float[InvSense.size()];
        y_QxData = new float[InvSense.size()];
        z_QxData = new float[InvSense.size()];
        //QxR*Qxs^-1;
        for(int i=0; i<InvSense.size();i++){
            QxRxInvSense1.add(quatMul(QxR.get(i), InvSense.get(i)));
            QxRxInvSense.add(QxRxInvSense1.get(i).clone());
            /////////////////////////
            w_QxData[i] = QxRxInvSense.get(i)[0];
            x_QxData[i] = QxRxInvSense.get(i)[1];
            y_QxData[i] = QxRxInvSense.get(i)[2];
            z_QxData[i] = QxRxInvSense.get(i)[3];
            /////////////////////////
        }
        //log QxData
        w_LogQxData = new float[QxRxInvSense.size()];
        x_LogQxData = new float[QxRxInvSense.size()];
        y_LogQxData = new float[QxRxInvSense.size()];
        z_LogQxData = new float[QxRxInvSense.size()];
        for(int i=0; i<QxRxInvSense.size(); i++){
            QxRxInvSenseLog1.add(logMap(QxRxInvSense.get(i)));
            QxRxInvSenseLog.add(QxRxInvSenseLog1.get(i).clone());
            w_LogQxData[i]=0;
            x_LogQxData[i]=QxRxInvSenseLog.get(i)[0];
            y_LogQxData[i]=QxRxInvSenseLog.get(i)[1];
            z_LogQxData[i]=QxRxInvSenseLog.get(i)[2];
        }


        // Dlugosc wektora i calka
        for(int i =0; i<x_LogQxData.length;i++){
            V_Length.add((float)Math.sqrt((double)(x_LogQxData[i]*x_LogQxData[i]
                    +y_LogQxData[i]*y_LogQxData[i]+z_LogQxData[i]*z_LogQxData[i])));
            error+=V_Length.get(i);
        }

        int min;
        int max;
        float errorMax;
        float errorMin;
        min = V_Length.indexOf(Collections.min((V_Length)));
        max = V_Length.indexOf(Collections.max((V_Length)));
        errorMin = V_Length.get(min);
        errorMax = V_Length.get(max);
        min =0;


        for (int i=0; i<V_Length.size(); i++){
            min = V_Length.indexOf(Collections.min((V_Length)));
            V_Length_sort.add(V_Length.get(min));
            V_Length.remove(min);
        }
        min =0;

    }


    public void buttonOnClick(View v){
        Button button = (Button) v;
        ((Button) v).setText("Iteracja: "+Float.toString(licznik_przycisk));

        //wyswietlanie orient
        if(t_orient_int.contains(licznik_przycisk)) {
            TextView Magnetx = (TextView)findViewById(R.id.Orientx);
            TextView Magnety = (TextView)findViewById(R.id.Orienty);
            TextView Magnetz = (TextView)findViewById(R.id.Orientz);
            Magnetx.setText(Float.toString(x_orient.get(t_orient_int.indexOf(licznik_przycisk))));
            Magnety.setText(Float.toString(y_orient.get(t_orient_int.indexOf(licznik_przycisk))));
            Magnetz.setText(Float.toString(z_orient.get(t_orient_int.indexOf(licznik_przycisk))));

            //obliczanie pitch roll yaw dla orient
            Pitch_orient.add((float) Math.atan2(y_orient.get(t_orient_int.indexOf(licznik_przycisk)),
                    z_orient.get(t_orient_int.indexOf(licznik_przycisk))) * 180 / (float) PI);
            Roll_orient.add((float) Math.atan2(x_orient.get(t_orient_int.indexOf(licznik_przycisk)),
                    z_orient.get(t_orient_int.indexOf(licznik_przycisk))) * 180 / (float) PI);
            Yaw_orient.add((float) Math.atan2(x_orient.get(t_orient_int.indexOf(licznik_przycisk)),
                    y_orient.get(t_orient_int.indexOf(licznik_przycisk))) * 180 / (float) PI);
        }
        else{
            if(!Pitch_orient.isEmpty()){
                Pitch_orient.add(Pitch_orient.get(Pitch_orient.size()-1));
                Roll_orient.add(Roll_orient.get(Roll_orient.size()-1));
                Yaw_orient.add(Yaw_orient.get(Yaw_orient.size()-1));
            }

        }


        //wyswietlanie acc
        if(t_acc_int.contains(licznik_przycisk)){
            TextView myTextView = (TextView)findViewById(R.id.myTextView3);
            TextView myTextView2 = (TextView)findViewById(R.id.myTextView4);
            TextView myTextView3 = (TextView)findViewById(R.id.myTextView5);
            myTextView.setText(Float.toString(x_acc.get(t_acc_int.indexOf(licznik_przycisk))));
            myTextView2.setText(Float.toString(y_acc.get(t_acc_int.indexOf(licznik_przycisk))));
            myTextView3.setText(Float.toString(z_acc.get(t_acc_int.indexOf(licznik_przycisk))));
        }


        //wyswietlanie gyro
        if(t_gyro_int.contains(licznik_przycisk)) {
            TextView Gyrox = (TextView)findViewById(R.id.Gyrox);
            TextView Gyroy = (TextView)findViewById(R.id.Gyroy);
            TextView Gyroz = (TextView)findViewById(R.id.Gyroz);
            Gyrox.setText(Float.toString(x_gyro.get(t_gyro_int.indexOf(licznik_przycisk))));
            Gyroy.setText(Float.toString(y_gyro.get(t_gyro_int.indexOf(licznik_przycisk))));
            Gyroz.setText(Float.toString(z_gyro.get(t_gyro_int.indexOf(licznik_przycisk))));

            //obliczanie pitch roll yaw dla gyro
            gyro_sum_x+=x_gyro.get(t_gyro_int.indexOf(licznik_przycisk));
            gyro_sum_y+=y_gyro.get(t_gyro_int.indexOf(licznik_przycisk));
            gyro_sum_z+=z_gyro.get(t_gyro_int.indexOf(licznik_przycisk));

            Pitch_gyro.add(gyro_sum_x);
            Roll_gyro.add(gyro_sum_y);
            Yaw_gyro.add(gyro_sum_z);
        }
        else{
            if (!Pitch_gyro.isEmpty()){
                Pitch_gyro.add(Pitch_gyro.get(Pitch_gyro.size()-1));
                Roll_gyro.add(Roll_gyro.get(Roll_gyro.size()-1));
                Yaw_gyro.add(Yaw_gyro.get(Yaw_gyro.size()-1));
            }
        }



        //obliczanie pitch roll yaw
        if (Pitch_gyro.isEmpty()==false && Pitch_orient.isEmpty()==false) {
            Pitch.add((float) (0.98 * Pitch_gyro.get(Pitch_gyro.size() - 1)
                    + 0.02 * Pitch_orient.get(Pitch_orient.size() - 1)));
            Roll.add((float) (0.98 * Roll_gyro.get(Roll_gyro.size() - 1)
                    + 0.02 * Roll_orient.get(Roll_orient.size() - 1)));
            Yaw.add((float) (0.98 * Yaw_gyro.get(Yaw_gyro.size() - 1)
                    + 0.02 * Yaw_orient.get(Yaw_orient.size() - 1)));
        }

        //wyswietlanie pitch roll yaw
        TextView Pitch_txt = (TextView)findViewById(R.id.Pitch);
        TextView Roll_txt = (TextView)findViewById(R.id.Roll);
        TextView Yaw_txt = (TextView)findViewById(R.id.Yaw);
        if (Pitch.isEmpty()==false) {
            Pitch_txt.setText(Float.toString(Pitch.get(Pitch.size() - 1)));
            Roll_txt.setText(Float.toString(Roll.get(Roll.size() - 1)));
            Yaw_txt.setText(Float.toString(Yaw.get(Yaw.size() - 1)));
        }

        //wyswietlanie quaternion gyro
        if(t_gyro_int.contains(licznik_przycisk)) {
            TextView QuatGyroWtxt = (TextView)findViewById(R.id.QuatGyroW);
            TextView QuatGyroXtxt = (TextView)findViewById(R.id.QuatGyroX);
            TextView QuatGyroYtxt = (TextView)findViewById(R.id.QuatGyroY);
            TextView QuatGyroZtxt = (TextView)findViewById(R.id.QuatGyroZ);
            QuatGyroWtxt.setText(Float.toString(QuatGyroW[t_gyro_int.indexOf(licznik_przycisk)]));
            QuatGyroXtxt.setText(Float.toString(QuatGyroX[t_gyro_int.indexOf(licznik_przycisk)]));
            QuatGyroYtxt.setText(Float.toString(QuatGyroY[t_gyro_int.indexOf(licznik_przycisk)]));
            QuatGyroZtxt.setText(Float.toString(QuatGyroZ[t_gyro_int.indexOf(licznik_przycisk)]));


            //wyswietlanie log quaternion gyro
            TextView LogGyroWtxt = (TextView)findViewById(R.id.LogGyroW);
            TextView LogGyroXtxt = (TextView)findViewById(R.id.LogGyroX);
            TextView LogGyroYtxt = (TextView)findViewById(R.id.LogGyroY);
            TextView LogGyroZtxt = (TextView)findViewById(R.id.LogGyroZ);
            LogGyroWtxt.setText(Float.toString(w[t_gyro_int.indexOf(licznik_przycisk)]));
            LogGyroXtxt.setText(Float.toString(x[t_gyro_int.indexOf(licznik_przycisk)]));
            LogGyroYtxt.setText(Float.toString(y[t_gyro_int.indexOf(licznik_przycisk)]));
            LogGyroZtxt.setText(Float.toString(z[t_gyro_int.indexOf(licznik_przycisk)]));

        }

        //wyswietlanie quaternion orient
        if(t_orient_int.contains(licznik_przycisk)){
            TextView QuatOrientWtxt = (TextView)findViewById(R.id.QuatOrientW);
            TextView QuatOrientXtxt = (TextView)findViewById(R.id.QuatOrientX);
            TextView QuatOrientYtxt = (TextView)findViewById(R.id.QuatOrientY);
            TextView QuatOrientZtxt = (TextView)findViewById(R.id.QuatOrientZ);
            QuatOrientWtxt.setText(Float.toString(quaternion_acm.get(t_orient_int.indexOf(licznik_przycisk))[0]));
            QuatOrientXtxt.setText(Float.toString(quaternion_acm.get(t_orient_int.indexOf(licznik_przycisk))[1]));
            QuatOrientYtxt.setText(Float.toString(quaternion_acm.get(t_orient_int.indexOf(licznik_przycisk))[2]));
            QuatOrientZtxt.setText(Float.toString(quaternion_acm.get(t_orient_int.indexOf(licznik_przycisk))[3]));


            //wyswietlanie log quaternion orient
            TextView LogOrientWtxt = (TextView)findViewById(R.id.LogOrientW);
            TextView LogOrientXtxt = (TextView)findViewById(R.id.LogOrientX);
            TextView LogOrientYtxt = (TextView)findViewById(R.id.LogOrientY);
            TextView LogOrientZtxt = (TextView)findViewById(R.id.LogOrientZ);
            LogOrientWtxt.setText(Float.toString(a[t_orient_int.indexOf(licznik_przycisk)]));
            LogOrientXtxt.setText(Float.toString(b[t_orient_int.indexOf(licznik_przycisk)]));
            LogOrientYtxt.setText(Float.toString(c[t_orient_int.indexOf(licznik_przycisk)]));
            LogOrientZtxt.setText(Float.toString(d[t_orient_int.indexOf(licznik_przycisk)]));
        }


        //wyswietlanie wyjscia filtru

        TextView QuatFilterW = (TextView)findViewById(R.id.W_QFilter);
        TextView QuatFilterX = (TextView)findViewById(R.id.X_QFilter);
        TextView QuatFilterY = (TextView)findViewById(R.id.Y_QFilter);
        TextView QuatFilterZ = (TextView)findViewById(R.id.Z_QFilter);
        QuatFilterW.setText(Float.toString(w_filter[licznik_przycisk]));
        QuatFilterX.setText(Float.toString(x_filter[licznik_przycisk]));
        QuatFilterY.setText(Float.toString(y_filter[licznik_przycisk]));
        QuatFilterZ.setText(Float.toString(z_filter[licznik_przycisk]));

        //wyswietlanie exp quaternion orient
        TextView QuatFilterExpW = (TextView)findViewById(R.id.W_Exp);
        TextView QuatFilterExpX = (TextView)findViewById(R.id.X_Exp);
        TextView QuatFilterExpY = (TextView)findViewById(R.id.Y_Exp);
        TextView QuatFilterExpZ = (TextView)findViewById(R.id.Z_Exp);
        QuatFilterExpW.setText(Float.toString(w_filter_exp[licznik_przycisk]));
        QuatFilterExpX.setText(Float.toString(x_filter_exp[licznik_przycisk]));
        QuatFilterExpY.setText(Float.toString(y_filter_exp[licznik_przycisk]));
        QuatFilterExpZ.setText(Float.toString(z_filter_exp[licznik_przycisk]));





        licznik_przycisk++;

    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        // Inflate the menu; this adds items to the action bar if it is present.
        getMenuInflater().inflate(R.menu.menu_main, menu);

        return true;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        // Handle action bar item clicks here. The action bar will
        // automatically handle clicks on the Home/Up button, so long
        // as you specify a parent activity in AndroidManifest.xml.
        int id = item.getItemId();

        //noinspection SimplifiableIfStatement
        if (id == R.id.action_settings) {
            return true;
        }

        return super.onOptionsItemSelected(item);
    }
}
