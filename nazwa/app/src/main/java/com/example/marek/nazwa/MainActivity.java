package com.example.marek.nazwa;

import android.hardware.SensorManager;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.support.v7.widget.Toolbar;
import android.view.View;
import android.view.Menu;
import android.view.MenuItem;
import android.widget.Button;
import android.widget.TextView;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class MainActivity extends AppCompatActivity {
    float[] acc;
    float[] data;
    String[] acc_string;
    String[] data_string;
    float[] gyro;
    float[] magnet;
    float[] orient;
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
    ArrayList<Float> t_magnet = new ArrayList<Float>();
    ArrayList<Float> x_magnet = new ArrayList<Float>();
    ArrayList<Float> y_magnet = new ArrayList<Float>();
    ArrayList<Float> z_magnet = new ArrayList<Float>();
    ArrayList<Float> Pitch_acc = new ArrayList<Float>();
    ArrayList<Float> Roll_acc = new ArrayList<Float>();
    ArrayList<Float> Yaw_acc = new ArrayList<Float>();
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
    ArrayList<float[]> Qacm = new ArrayList<>();
    ArrayList<float[]> v = new ArrayList<>();
    ArrayList<Float> t_data = new ArrayList<Float>();
    ArrayList<Float> w_data = new ArrayList<Float>();
    ArrayList<Float> x_data = new ArrayList<Float>();
    ArrayList<Float> y_data = new ArrayList<Float>();
    ArrayList<Float> z_data = new ArrayList<Float>();
    ArrayList<Integer> t_data_int = new ArrayList<Integer>();
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
    ArrayList<float[]> quaternion_gyro1 = new ArrayList<>();

    float[] OrientQuat =null;
    float[] qacm = new float [3];
    float[] qgyro = new float [3];
    float[] GyroQuat =new float [4];
    float[] Qg0=null;
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

    int licznik_przycisk=0;
    Float gyro_sum_x= Float.valueOf(0);
    Float gyro_sum_y=Float.valueOf(0);
    Float gyro_sum_z=Float.valueOf(0);
    static final double ms2s = 1e-3;
    double PI = 3.14159265359;
    private static final float eps = 1e-9f;




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



        try {
            reader = new BufferedReader(
                    new InputStreamReader(getAssets().open("acc.log")));

            reader1 = new BufferedReader(
                    new InputStreamReader(getAssets().open("gyro.log")));
            reader2 = new BufferedReader(
                    new InputStreamReader(getAssets().open("mag.log")));

            reader3 = new BufferedReader(
                    new InputStreamReader(getAssets().open("orient.log")));
            reader4 = new BufferedReader(
                    new InputStreamReader(getAssets().open("Xsense.data")));

            // do reading, usually loop until end of file reading

            while ((mLine = reader.readLine()) != null) {
                //process line
                acc_string = mLine.split(" ");
                acc = new float[acc_string.length];

                for (int i=0; i<acc_string.length; i++){
                    acc[i]=Float.valueOf(acc_string[i]);
                }
                t_acc.add(acc[0]);
                t_acc_int.add((int)acc[0]);
                x_acc.add(acc[1]);
                y_acc.add(acc[2]);
                z_acc.add(acc[3]);
            }
            while ((mLine1 = reader1.readLine()) != null) {
                //process line
                String[] gyro_string = mLine1.split(" ");
                gyro = new float[gyro_string.length];

                for (int i=0; i<gyro_string.length; i++){
                    gyro[i]=Float.valueOf(gyro_string[i]);
                }

                t_gyro.add(gyro[0]);
                t_gyro_int.add((int)gyro[0]);
                x_gyro.add(gyro[1]);
                y_gyro.add(gyro[2]);
                z_gyro.add(gyro[3]);

            }
            for (int i=0; i<t_gyro_int.size();i++){
                qgyro[0]=(float)x_gyro.get(i);
                qgyro[1]=(float)y_gyro.get(i);
                qgyro[2]=(float)z_gyro.get(i);

                float dt;
                if (t_gyro_int.size()==1) {
                    dt=1*(float)ms2s;
                }
                else{
                    dt = (t_gyro_int.get(t_gyro_int.size() - 1) - t_gyro_int.get(t_gyro_int.size() - 2))*(float)ms2s;
                }
                compQuatFromGyro(qgyro, dt , GyroQuat);
                if(quaternion_gyro1.isEmpty()){
                    quaternion_gyro1.add(GyroQuat);
                    quaternion_gyro.add(quaternion_gyro1.get(i).clone());
                }
                else {
                    quaternion_gyro1.add(quaternion_gyro1.size() - 1, GyroQuat);
                    quaternion_gyro.add(quaternion_gyro1.get(i).clone());
                }
            }
            w = new float[quaternion_gyro.size()];
            x = new float[quaternion_gyro.size()];
            y = new float[quaternion_gyro.size()];
            z = new float[quaternion_gyro.size()];
            QuatGyroW = new float[quaternion_gyro.size()];
            QuatGyroX = new float[quaternion_gyro.size()];
            QuatGyroY = new float[quaternion_gyro.size()];
            QuatGyroZ = new float[quaternion_gyro.size()];
            for(int i=0;i<quaternion_gyro.size();i++) {
                float[] pomocnicza = null;
                pomocnicza = quaternion_gyro.get(i);
                w[i] = pomocnicza[0];
                x[i] = pomocnicza[1];
                y[i] = pomocnicza[2];
                z[i] = pomocnicza[3];

                if (i==0){
                    Qg0=quaternion_gyro.get(i);
                }
                float[]pomocnicza1=null;
                if(i>0){
                    pomocnicza1=quaternion_gyro.get(i-1);
                    w[i]=pomocnicza1[0]*pomocnicza[0]-(pomocnicza1[1]*pomocnicza[1]+pomocnicza1[2]*pomocnicza[2]+pomocnicza1[3]*pomocnicza[3]);

                    x[i]=pomocnicza1[0]*pomocnicza[1]+pomocnicza[0]*pomocnicza1[1]+pomocnicza1[2]*pomocnicza[3]-pomocnicza1[3]*pomocnicza[2];
                    y[i]=pomocnicza1[0]*pomocnicza[2]+pomocnicza[0]*pomocnicza1[2]+pomocnicza1[3]*pomocnicza[1]-pomocnicza1[1]*pomocnicza[3];
                    z[i]=pomocnicza1[0]*pomocnicza[3]+pomocnicza[0]*pomocnicza1[3]+pomocnicza1[1]*pomocnicza[2]-pomocnicza1[2]*pomocnicza[1];
                    QuatGyroW[i] = w[i];
                    QuatGyroX[i] = x[i];
                    QuatGyroY[i] = y[i];
                    QuatGyroZ[i] = z[i];

                }
                double Length;

                Length=Math.sqrt((double)(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]));
                Length=Math.atan(Length/w[i]);
                w[i]=0;
                x[i]*=Length;
                y[i]*=Length;
                z[i]*=Length;
            }




            while ((mLine2 = reader2.readLine()) != null) {
                //process line
                String[] magnet_string = mLine2.split(" ");
                magnet=new float[magnet_string.length];

                for (int i=0; i<magnet_string.length; i++){
                    magnet[i]=Float.valueOf(magnet_string[i]);
                }
                t_magnet.add(magnet[0]);
                x_magnet.add(magnet[1]);
                y_magnet.add(magnet[2]);
                z_magnet.add(magnet[3]);
            }
            while ((mLine3 = reader3.readLine()) != null) {
                //process line
                String[] orient_string = mLine3.split(" ");
                orient=new float[orient_string.length];

                for (int i=0; i<orient_string.length; i++){
                    orient[i]=Float.valueOf(orient_string[i]);

                }


                t_orient.add(orient[0]);
                qacm[0]=(float)orient[1];
                qacm[1]=(float)orient[2];
                qacm[2]=(float)orient[3];
                t_orient_int.add((int)orient[0]);
                x_orient.add(orient[1]);
                y_orient.add(orient[2]);
                z_orient.add(orient[3]);
                if(OrientQuat==null)
                {
                    OrientQuat = new float[4];
                }
                //if(t_orient_int.contains(k))
                //{
                    SensorManager.getQuaternionFromVector(OrientQuat,qacm);
                    quaternion_acm1.add(OrientQuat);
                    quaternion_acm.add(quaternion_acm1.get(quaternion_acm1.size()-1).clone());
               // }
                //else{
                //    quaternion_acm1.add(quaternion_acm.get(quaternion_acm.size()-1));
                //    quaternion_acm.add(quaternion_acm1.get(quaternion_acm1.size()-1).clone());
               // }


                k++;
            }

            a = new float[quaternion_acm.size()];
            b = new float[quaternion_acm.size()];
            c = new float[quaternion_acm.size()];
            d = new float[quaternion_acm.size()];
            float [] wektorV=new float[3];
            for(int i=0;i<quaternion_acm.size();i++){
                float [] pomocnicza=null;
                pomocnicza=quaternion_acm.get(i);
                a[i]=pomocnicza[0];
                wektorV[0]=pomocnicza[1];
                wektorV[1]=pomocnicza[2];
                wektorV[2]=pomocnicza[3];
                b[i]=pomocnicza[1];
                c[i]=pomocnicza[2];
                d[i]=pomocnicza[3];
                v.add(wektorV);
                double Length;

                Length=Math.sqrt((double)(b[i]*b[i]+c[i]*c[i]+d[i]*d[i]));
                Length=Math.atan(Length/a[i]);
                a[i]=0;
                b[i]*=Length;
                c[i]*=Length;
                d[i]*=Length;
            }

            float waga=0.98f;
            w_filter = new float[a.length];
            x_filter = new float[a.length];
            y_filter = new float[a.length];
            z_filter = new float[a.length];
            w_filter_exp = new float[a.length];
            x_filter_exp = new float[a.length];
            y_filter_exp = new float[a.length];
            z_filter_exp = new float[a.length];


            for (int i=0; i<a.length ;i++){
                ///mnozenie z waga
                w_filter[i]=waga*w[i]+(1-waga)*a[i];
                x_filter[i]=waga*x[i]+(1-waga)*b[i];
                y_filter[i]=waga*y[i]+(1-waga)*c[i];
                z_filter[i]=waga*z[i]+(1-waga)*d[i];
                w_filter_exp[i] = waga*w[i]+(1-waga)*a[i];
                x_filter_exp[i] = waga*x[i]+(1-waga)*b[i];
                y_filter_exp[i] = waga*y[i]+(1-waga)*c[i];
                z_filter_exp[i] = waga*z[i]+(1-waga)*d[i];

                ///exp
                double Mul;
                double Length;
                Length=Math.sqrt((double)(x_filter_exp[i]*x_filter_exp[i]+y_filter_exp[i]*y_filter_exp[i]+z_filter_exp[i]*z_filter_exp[i]));
                if (Length > 1.0e-4)
                    Mul = Math.sin(Length)/Length;
                else
                    Mul = 1.0;
                w_filter_exp[i] = (float)Math.cos(Length);
                x_filter_exp[i]*=Mul;
                y_filter_exp[i]*=Mul;
                z_filter_exp[i]*=Mul;

            }

            while ((mLine4 = reader4.readLine()) != null) {
                //process line
                data_string = mLine4.split(" ");
                data = new float[data_string.length];

                for (int i=0; i<data_string.length; i++){
                    data[i]=Float.valueOf(data_string[i]);
                }
                t_data.add(data[0]);
                //t_acc_int.add((int)acc[0]);
                w_data.add(data[4]);
                x_data.add(data[5]);
                y_data.add(data[6]);
                z_data.add(data[7]);
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

            //Q_filter * Rotation matrix Quaternion
            w_QxR = new float[w_filter_exp.length];
            x_QxR = new float[w_filter_exp.length];
            y_QxR = new float[w_filter_exp.length];
            z_QxR = new float[w_filter_exp.length];
            for(int i=0; i<w_filter_exp.length;i++){
                w_QxR[i]=w_filter_exp[i]*qw - x_filter_exp[i]*qx - y_filter_exp[i]*qy - z_filter_exp[i]*qz;
                x_QxR[i]=w_filter_exp[i]*qx + x_filter_exp[i]*qw + y_filter_exp[i]*qz - z_filter_exp[i]*qy;
                y_QxR[i]=w_filter_exp[i]*qy + y_filter_exp[i]*qw + z_filter_exp[i]*qx - x_filter_exp[i]*qz;
                z_QxR[i]=w_filter_exp[i]*qz + z_filter_exp[i]*qw + x_filter_exp[i]*qy - y_filter_exp[i]*qx;
            }


         //invert quaternion data (Xsense)
            for (int i=0; i<t_data.size();i++){
                x_data.set(i, -x_data.get(i));
                y_data.set(i, -y_data.get(i));
                z_data.set(i, -z_data.get(i));
            }

         //QxR*Qxs^-1;
            w_QxData = new float[w_QxR.length];
            x_QxData = new float[w_QxR.length];
            y_QxData = new float[w_QxR.length];
            z_QxData = new float[w_QxR.length];
            for(int i=0; i<w_data.size();i++){
                w_QxData[i]=w_QxR[i]*w_data.get(i) - x_QxR[i]*x_data.get(i) - y_QxR[i]*y_data.get(i) - z_QxR[i]*z_data.get(i);
                x_QxData[i]=w_QxR[i]*x_data.get(i) + x_QxR[i]*w_data.get(i) + y_QxR[i]*z_data.get(i) - z_QxR[i]*y_data.get(i);
                y_QxData[i]=w_QxR[i]*y_data.get(i) + y_QxR[i]*w_data.get(i) + z_QxR[i]*x_data.get(i) - x_QxR[i]*z_data.get(i);
                z_QxData[i]=w_QxR[i]*z_data.get(i) + z_QxR[i]*w_data.get(i) + x_QxR[i]*y_data.get(i) - y_QxR[i]*x_data.get(i);
            }

            //log QxData
           // w_LogQxData = new float[w_QxData.length];
           // x_LogQxData = new float[w_QxData.length];
           // y_LogQxData = new float[w_QxData.length];
           // z_LogQxData = new float[w_QxData.length];
           // for (int i=0; i<w_QxData.length;i++){
             //   double Length;

             //   Length=Math.sqrt((double)(x_QxData[i]*x_QxData[i]+y_QxData[i]*y_QxData[i]+z_QxData[i]*z_QxData[i]));
             //   Length=Math.atan(Length/w_QxData[i]);
             //   w_LogQxData[i]=0;
             //   x_LogQxData[i]*=Length;
             //   y_LogQxData[i]*=Length;
             //   z_LogQxData[i]*=Length;
            //}


            /*double LengthV;
            LengthV=Math.sqrt((double)(b[1]*b[1]+c[1]*c[1]+d[1]*d[1]));
            double LengthQ;
            LengthQ=Math.sqrt((double)(a[1]*a[1]+b[1]*b[1]+c[1]*c[1]+d[1]*d[1]));*/


        } catch (IOException e) {
            //log the exception
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    //log the exception
                }
            }
        }



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
        }

        //wyswietlanie acc
        if(t_acc_int.contains(licznik_przycisk)){
            TextView myTextView = (TextView)findViewById(R.id.myTextView3);
            TextView myTextView2 = (TextView)findViewById(R.id.myTextView4);
            TextView myTextView3 = (TextView)findViewById(R.id.myTextView5);
            myTextView.setText(Float.toString(x_acc.get(t_acc_int.indexOf(licznik_przycisk))));
            myTextView2.setText(Float.toString(y_acc.get(t_acc_int.indexOf(licznik_przycisk))));
            myTextView3.setText(Float.toString(z_acc.get(t_acc_int.indexOf(licznik_przycisk))));

            //obliczanie pitch roll yaw dla acc
            Pitch_acc.add((float) Math.atan2(y_acc.get(t_acc_int.indexOf(licznik_przycisk)), z_acc.get(t_acc_int.indexOf(licznik_przycisk))) * 180 / (float) PI);
            Roll_acc.add((float) Math.atan2(x_acc.get(t_acc_int.indexOf(licznik_przycisk)), z_acc.get(t_acc_int.indexOf(licznik_przycisk))) * 180 / (float) PI);
            Yaw_acc.add((float) Math.atan2(x_acc.get(t_acc_int.indexOf(licznik_przycisk)), y_acc.get(t_acc_int.indexOf(licznik_przycisk))) * 180 / (float) PI);
        }
        else{
            if(!Pitch_acc.isEmpty()){
                Pitch_acc.add(Pitch_acc.get(Pitch_acc.size()-1));
                Roll_acc.add(Roll_acc.get(Roll_acc.size()-1));
                Yaw_acc.add(Yaw_acc.get(Yaw_acc.size()-1));
            }

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
        if (Pitch_gyro.isEmpty()==false && Pitch_acc.isEmpty()==false) {
            Pitch.add((float) (0.98 * Pitch_gyro.get(Pitch_gyro.size() - 1) + 0.02 * Pitch_acc.get(Pitch_acc.size() - 1)));
            Roll.add((float) (0.98 * Roll_gyro.get(Roll_gyro.size() - 1) + 0.02 * Roll_acc.get(Roll_acc.size() - 1)));
            Yaw.add((float) (0.98 * Yaw_gyro.get(Yaw_gyro.size() - 1) + 0.02 * Yaw_acc.get(Yaw_acc.size() - 1)));
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
            //LogGyroWtxt.setText(Float.toString(w[t_gyro_int.indexOf(licznik_przycisk)]));
            //LogGyroXtxt.setText(Float.toString(x[t_gyro_int.indexOf(licznik_przycisk)]));
            //LogGyroYtxt.setText(Float.toString(y[t_gyro_int.indexOf(licznik_przycisk)]));
            //LogGyroZtxt.setText(Float.toString(z[t_gyro_int.indexOf(licznik_przycisk)]));
            LogGyroWtxt.setText(Float.toString(w[licznik_przycisk]));
            LogGyroXtxt.setText(Float.toString(x[licznik_przycisk]));
            LogGyroYtxt.setText(Float.toString(y[licznik_przycisk]));
            LogGyroZtxt.setText(Float.toString(z[licznik_przycisk]));

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
            //LogOrientWtxt.setText(Float.toString(a[t_orient_int.indexOf(licznik_przycisk)]));
            //LogOrientXtxt.setText(Float.toString(b[t_orient_int.indexOf(licznik_przycisk)]));
            //LogOrientYtxt.setText(Float.toString(c[t_orient_int.indexOf(licznik_przycisk)]));
            //LogOrientZtxt.setText(Float.toString(d[t_orient_int.indexOf(licznik_przycisk)]));
            LogOrientWtxt.setText(Float.toString(a[licznik_przycisk]));
            LogOrientXtxt.setText(Float.toString(b[licznik_przycisk]));
            LogOrientYtxt.setText(Float.toString(c[licznik_przycisk]));
            LogOrientZtxt.setText(Float.toString(d[licznik_przycisk]));
        }


        //wyswietlanie wyjscia filtru
        if ((t_orient_int.contains(licznik_przycisk)) || (t_gyro_int.contains(licznik_przycisk))){
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
        }




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
