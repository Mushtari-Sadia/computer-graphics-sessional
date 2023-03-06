#include <bits/stdc++.h> // header in local directory
#include <iostream> // header in standard library
#include <fstream>
#define pi (2*acos(0.0))
#include "bitmap_image.hpp"
using namespace std;

struct point
{
	double x,y,z;
};

struct color
{
    double red,green,blue;
};
struct triangle
{
    struct point points[3];
    struct color clr;
};

double **z_buffer;
struct color **frame_buffer;

struct point eye,look,up,fov;
double aspectRatio,near,far;
double **V,**P;
ofstream out("stage1.txt");
ofstream out2("stage2.txt");
ofstream out3("stage3.txt");
ofstream out4;

int triangle_count=0;
struct triangle *triangles;

double** gen_identity_matrix()
{
    double** matrix = new double*[4];
    for(int i = 0; i < 4; i++){
        matrix[i] = new double[4];
        for(int j = 0; j < 4; j++){
            if(i==j)
            {
                matrix[i][j]=1;
            }
            else
            {
                matrix[i][j]=0;
            }
        }
    }
    return matrix;
}


double** product(double **a, double **b)
{
    //cout << "inside product" << endl;
    double **mult = new double*[4];
    for(int i = 0; i < 4; i++){
            mult[i] = new double[4];
            for(int j = 0; j < 4; ++j)
            {
                mult[i][j]=0;
            }
        }

    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
            for(int k = 0; k < 4; ++k)
            {
                mult[i][j] += a[i][k] * b[k][j];
            }
    return mult;
}

struct point R(struct point v1,struct point v2,double angle)
{
    struct point ans;
    angle = -angle*pi/180;
    double dot = (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
    ans.x = cos(angle)*v1.x + (1-cos(angle))*dot*v2.x + sin(angle)*(v1.y * v2.z - v1.z * v2.y);
    ans.y = cos(angle)*v1.y + (1-cos(angle))*dot*v2.y + sin(angle)*(-1)*(v1.x * v2.z - v1.z * v2.x);
    ans.z = cos(angle)*v1.z + (1-cos(angle))*dot*v2.z + sin(angle)*(v1.x* v2.y - v1.y * v2.x);
    return ans;
}

struct point normalize(struct point a)
{
    struct point n;
    double length = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    n.x = a.x/length;
    n.y = a.y/length;
    n.z = a.z/length;
    return n;
};


double** generate_rotation_matrix(struct point oa,double angle)
{
    struct point a= normalize(oa);
    struct point c1,c2,c3;
    c1 = R((struct point){1,0,0},a,angle);
    c2 = R((struct point){0,1,0},a,angle);
    c3 = R((struct point){0,0,1},a,angle);

    double** matrix = gen_identity_matrix();

    matrix[0][0] = c1.x;
    matrix[0][1] = c2.x;
    matrix[0][2] = c3.x;

    matrix[1][0] = c1.y;
    matrix[1][1] = c2.y;
    matrix[1][2] = c3.y;

    matrix[2][0] = c1.z;
    matrix[2][1] = c2.z;
    matrix[2][2] = c3.z;

    return matrix;

}


struct point TransformPoint(double **matrix,struct point P)
{
    struct point transformedPoint;
    //transformation
    transformedPoint.x = matrix[0][0]*P.x + matrix[0][1]*P.y + matrix[0][2]*P.z + matrix[0][3]*1 ;
    transformedPoint.y = matrix[1][0]*P.x + matrix[1][1]*P.y + matrix[1][2]*P.z + matrix[1][3]*1 ;
    transformedPoint.z = matrix[2][0]*P.x + matrix[2][1]*P.y + matrix[2][2]*P.z + matrix[2][3]*1 ;
    double w = matrix[3][0]*P.x + matrix[3][1]*P.y + matrix[3][2]*P.z + matrix[3][3]*1 ;
    if(w!=1)
    {
        transformedPoint.x /=w;
        transformedPoint.y /=w;
        transformedPoint.z /=w;
    }
    return transformedPoint;
};

void modeling_transformation()
{
    string line;


    stack<double**> stk;
    stack<stack<double**>> states;
    double** matrix = gen_identity_matrix();

    stk.push(matrix);

    while(std::getline(std::cin, line))
    {
         if(line=="triangle")
         {
             for(int i=0;i<3;i++)
                {
                    struct point p;
                    cin >> line;
                    //cout << line << " ";
                    p.x = stod(line);

                    cin >> line;
                    //cout << line << " ";
                    p.y = stod(line);

                    cin >> line;
                    //cout << line << " ";
                    p.z = stod(line);

                    struct point TP = TransformPoint(stk.top(),p);
                    out << std::fixed;
                    out << std::setprecision(7);
                    out << TP.x << " " << TP.y << " " << TP.z << " " << endl;


                    struct point VTP = TransformPoint(V,TP);
                    out2 << std::fixed;
                    out2 << std::setprecision(7);
                    out2 << VTP.x << " " << VTP.y << " " << VTP.z << " " << endl;

                    struct point PTP = TransformPoint(P,VTP);
                    out3 << std::fixed;
                    out3 << std::setprecision(7);
                    out3 << PTP.x << " " << PTP.y << " " << PTP.z << " " << endl;


                }
                triangle_count++;
                out << endl;
                out2 << endl;
                out3 << endl;
         }

         else if(line=="translate")
         {
            double tx,ty,tz;
            cin >> tx >> ty >> tz;
            //cout << tx << ty << tz;
            //generate translation matrix
            double** tr = gen_identity_matrix();
            tr[0][3] = tx;
            tr[1][3] = ty;
            tr[2][3] = tz;
            //cout << "calling product" << endl;
            stk.push(product(stk.top(),tr));
            for(int k=0; k<4; k++)
            {
                delete [] tr[k];
            }
            delete [] tr;
         }

         else if(line=="scale")
         {
            double sx,sy,sz;
            cin >> sx >> sy >> sz;
            //cout << sx << sy << sz;
            //generate scale matrix
            double** tr = gen_identity_matrix();
            tr[0][0] = sx;
            tr[1][1] = sy;
            tr[2][2] = sz;
            //cout << "calling product" << endl;
            stk.push(product(stk.top(),tr));
            //cout << "pushed to stack" << endl;
            delete [] tr;
         }

         else if(line=="rotate")
         {
             double angle;
             struct point a;
             cin >> angle >> a.x >> a.y >> a.z;
             double **tr=generate_rotation_matrix(a,angle);
             stk.push(product(stk.top(),tr));
             delete [] tr;
         }
         else if(line=="push")
         {
            states.push(stk);
         }
         else if(line=="pop")
         {
             stk = states.top();
             states.pop();
         }
         else if(line=="end")
         {
             break;
         }
    }

}

struct point cross(struct point a, struct point b)
{
    struct point res;
    res.x = a.y*b.z - a.z*b.y;
    res.y = (-1)*(a.x*b.z - a.z*b.x);
    res.z = a.x*b.y - a.y*b.x;
    return res;
};

void view_transformation()
{
    struct point l,r,u;
    l.x = look.x - eye.x;
    l.y = look.y - eye.y;
    l.z = look.z - eye.z;

    l = normalize(l);
    r = cross(l,up);
    r = normalize(r);
    u = cross(r,l);

    double** T = gen_identity_matrix();
    T[0][3] = -eye.x;
    T[1][3] = -eye.y;
    T[2][3] = -eye.z;

    double** R = gen_identity_matrix();
    R[0][0] = r.x;
    R[0][1] = r.y;
    R[0][2] = r.z;

    R[1][0] = u.x;
    R[1][1] = u.y;
    R[1][2] = u.z;

    R[2][0] = -l.x;
    R[2][1] = -l.y;
    R[2][2] = -l.z;

    V=product(R,T);
    delete [] R;
    delete [] T;
}

void projection_transformation()
{
    //cout << "near " << near << endl;
    //cout << "far " << far << endl;
    //cout << "aspectRatio " << aspectRatio << endl;
    fov.x = fov.y*aspectRatio;
    double t = near*tan(fov.y/2*pi/180);
    double r = near*tan(fov.x/2*pi/180);

    P = gen_identity_matrix();
    P[0][0] = near/r;
    P[1][1] = near/t;
    P[2][2] = -(far+near)/(far-near);
    P[3][3] = 0;

    P[2][3] = -(2*far*near)/(far-near);
    P[3][2] = -1;
}

void printMatrix(double** matrix)
{
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            //cout << matrix[i][j] << " ";
        }
        //cout << endl;
    }
}


void find_intersection(double &intersect_1_x,double &intersect_2_x,struct triangle T, double scanline,double top_y,double dy)
{
    int count=0;
    for(int i=0;i<3;i++)
    {
        int point_a,point_b;
        if(i==0)
        {
            point_a=1;
            point_b=2;
        }
        else if(i==1)
        {
            point_a=0;
            point_b=2;
        }
        else
        {
            point_a=0;
            point_b=1;
        }

        //cout << "point_a" << point_a << endl;
        //cout << "point_b" << point_b << endl;
        struct point A = T.points[point_a];
        struct point B = T.points[point_b];

        double big_y = B.y;
        double small_y = A.y;
        if(A.y>B.y)
        {
            big_y = A.y;
            small_y = B.y;
        }

        ////cout << "big_y = " << big_y << endl;
        ////cout << "small_y = " << small_y << endl;

        big_y = round((top_y-big_y)/dy);
        small_y = round((top_y-small_y)/dy);

        ////cout << "big_y_row = " << big_y << endl;
        ////cout << "small_y_row = " << small_y << endl;

        if(scanline<=small_y && scanline>=big_y)
        {
            double a1,b1,a2,b2,c1,c2,determinant,x,y;
            a1 = B.y - A.y;
            b1 = A.x - B.x;
            c1 = a1*(A.x) + b1*(A.y);

            a2 = 0;
            b2 = 1;
            c2 = top_y - scanline*dy;

        //    //cout << "c2 = " << c2 << endl;

            determinant = a1*b2 - a2*b1;

            if (determinant == 0)
            {
                if(c1==c2)
                {
                    intersect_1_x = A.x;
                    intersect_2_x = B.x;
                    return;
                }
                else
                {
                    continue;
                }
            }

            else
            {
                x = (b2*c1 - b1*c2)/determinant;
                y = (a1*c2 - a2*c1)/determinant;
                if(count==0)
                {
                    intersect_1_x = x;
                }
                else
                {
                    intersect_2_x = x;
                    return;
                }
                 count++;
            }
        }
        else
            continue;
    }



};



void clip_scan_conversion()
{
    int screen_width,screen_height;
    double limit_x,limit_y,front_z,rear_z;

    ifstream in("config.txt");
    in >> screen_width >> screen_height;
    in >> limit_x;
    in >> limit_y;
    in >> front_z >> rear_z;

    //cout << "screen width " << screen_width << endl;
    //cout << "screen height " << screen_height << endl;
    //cout << "limit_x " << limit_x << endl;
    //cout << "limit_y " << limit_y << endl;
    //cout << "front_z " << front_z << endl;
    //cout << "rear_z " << rear_z << endl;


    ifstream in2("stage3.txt");
    triangles = new struct triangle[triangle_count];
    for(int i=0;i<triangle_count;i++)
    {
        for(int j=0;j<3;j++)
        {
            in2 >> triangles[i].points[j].x >> triangles[i].points[j].y >> triangles[i].points[j].z;
        }
        triangles[i].clr.red = 0 + rand() % (( 255 + 1 ) - 0);
        triangles[i].clr.green = 0 + rand() % (( 255 + 1 ) - 0);
        triangles[i].clr.blue = 0 + rand() % (( 255 + 1 ) - 0);
    }
    //cout << "triangles read" << endl;
    for(int i=0;i<triangle_count;i++)
    {
        for(int j=0;j<3;j++)
        {
            //cout << triangles[i].points[j].x <<" "<< triangles[i].points[j].y <<" "<< triangles[i].points[j].z << endl;

        }
        //cout << "color " << triangles[i].clr.red << " " << triangles[i].clr.green << " " <<triangles[i].clr.blue << endl;
    }


    //cout << "end" << endl;
    z_buffer = new double*[screen_width];
    for(int i = 0; i < screen_width; i++){
        z_buffer[i] = new double[screen_height];
        for(int j = 0; j < screen_height; j++){
                z_buffer[i][j] = rear_z;
        }
    }
    frame_buffer = new struct color*[screen_width];
    for(int i = 0; i < screen_width; i++){
        frame_buffer[i] = new struct color[screen_height];
        for(int j = 0; j < screen_height; j++){
                frame_buffer[i][j].red=0;
                frame_buffer[i][j].green=0;
                frame_buffer[i][j].blue=0;
        }
    }
//    //cout << "zbuffer" << endl;
//    for(int i = 0; i < screen_width; i++){
//        for(int j = 0; j < screen_height; j++){
//                //cout << z_buffer[i][j] << "\t";
//        }
//    }
//    //cout << endl << "frame buffer" << endl;
//    for(int i = 0; i < screen_width; i++){
//        for(int j = 0; j < screen_height; j++){
//               //cout << frame_buffer[i][j].red << "\t";
//        }
//    }
//    //cout << endl;

    double dx = (-limit_x-limit_x)/screen_width;
    double dy = (-limit_y-limit_y)/screen_height;
    double top_y = -limit_y-dy/2.0;
    double bottom_y = -top_y;
    double left_x = limit_x+dx/2.0;
    double right_x = -left_x;
    double max_y,min_y;
    int top_point,left_point,right_point;

    //cout << "dx= " << dx << endl;
    //cout << "dy= " << dy << endl;
    //cout << "top_y= " << top_y << endl;
    //cout << "bottom_y= " << bottom_y << endl;
    //cout << "left_x= " << left_x << endl;
    //cout << "right_x= " << right_x << endl;

    //cout << "triangle count " << triangle_count << endl;
    for(int i=0;i<triangle_count;i++)
    {
        max_y = triangles[i].points[0].y;
        min_y = triangles[i].points[0].y;

        top_point = 0;
        for(int j=0;j<3;j++)
        {
            if(triangles[i].points[j].y>max_y)
            {
                max_y = triangles[i].points[j].y;
                top_point=j;
            }
            if(triangles[i].points[j].y<min_y)
            {
                min_y = triangles[i].points[j].y;
            }
        }

        double top_scanline,bottom_scanline;

        if(max_y>=top_y)
        {
            top_scanline=0;
        }
        else
        {
            top_scanline = round((top_y-max_y)/dy);
        }

        if(min_y<=bottom_y)
        {
            bottom_scanline=screen_height-1;
        }
        else
        {
            bottom_scanline = screen_height-round((min_y-bottom_y)/dy)-1;
        }
        //cout << "max_y " << max_y << endl;
        //cout << "min_y " << min_y << endl;
        //cout << "top_scanline " << top_scanline << endl;
        //cout << "bottom_scanline " << bottom_scanline << endl;

        for(int row_no=(int)top_scanline;row_no<=(int)bottom_scanline;row_no++)
        {
            double left_int_col,right_int_col;
            double intersect_1_x,intersect_2_x;
            find_intersection(intersect_1_x,intersect_2_x,triangles[i],row_no,top_y,dy);
            double left_int_x = intersect_1_x;
            double right_int_x= intersect_2_x;
            if(intersect_1_x > intersect_2_x)
            {
                left_int_x = intersect_2_x;
                right_int_x= intersect_1_x;
            }

            //cout << "scanline " << row_no << endl;
            //cout << "left_int_x " << left_int_x << endl;
            //cout << "right_int_x " << right_int_x << endl;

            //clipping
            if(left_int_x<=left_x)
            {
//                left_int_x = round(left_x);
                  left_int_col = 0;
            }
            else
            {
                //find left intersecting col
                left_int_col = round((left_int_x-left_x)/dx);
            }

            //clipping
            if(right_int_x>=right_x)
            {
//                right_int_col = round(right_x);
                right_int_col = screen_width-1;
            }
            else
            {
                //find right intersecting col
                right_int_col = screen_width-1-round((right_x-right_int_x)/dx);
            }

            //cout << "left_int_col" << left_int_col << endl;
            //cout << "right_int_col" << right_int_col << endl;

            //find point1,point2,point3
            int sum,points[3]; //0 if top side of scanline,1 if bottom
            double ys = top_y - row_no*dy;
            sum=0;
            for(int k=0;k<3;k++)
            {
                if(triangles[i].points[k].y>=ys)
                    points[k]=0;
                else
                {
                    points[k]=1;
                    sum++;
                }
            }
            if(sum==2)
            {
                for(int k=0;k<3;k++)
                {
                    if(points[k]==0)
                        top_point=k;
                }
            }
            else if(sum==1)
            {
                for(int k=0;k<3;k++)
                {
                    if(points[k]==1)
                        top_point=k;
                }
            }

            //find z value
            double x1 = triangles[i].points[top_point].x;
            double y1 = triangles[i].points[top_point].y;
            double z1 = triangles[i].points[top_point].z;

            if(top_point==0)
            {
                left_point=1;
                right_point=2;
            }
            else if(top_point==1)
            {
                left_point=0;
                right_point=2;
            }
            else
            {
                left_point=0;
                right_point=1;
            }

            if(triangles[i].points[left_point].x>triangles[i].points[right_point].x)
            {
                int temp = left_point;
                left_point = right_point;
                right_point = temp;
            }

            double x2 = triangles[i].points[left_point].x;
            double y2 = triangles[i].points[left_point].y;
            double z2 = triangles[i].points[left_point].z;

            double x3 = triangles[i].points[right_point].x;
            double y3 = triangles[i].points[right_point].y;
            double z3 = triangles[i].points[right_point].z;

//            double ys= top_y-row_no*dy;
            double xa=left_int_x;
            double xb=right_int_x;
            double za,zb,zp;

            za = z1 + ((ys-y1)/(y2-y1))*(z2-z1);
            zb = z1 + ((ys-y1)/(y3-y1))*(z3-z1);

            double constant_term = (dx/(xb-xa))*(zb-za);

            for(int col_no=(int)left_int_col;col_no<=(int)right_int_col;col_no+=1)
            {

                if(col_no == (int)left_int_col)
                {
                    double xp=left_x+col_no*dx;
                    zp = za+((xp-xa)/(xb-xa))*(zb-za);

                }
                else
                    zp+=constant_term;




                if(zp>front_z && zp<z_buffer[row_no][col_no] && row_no>(int)top_scanline) {
                    z_buffer[row_no][col_no] = zp;
                    frame_buffer[row_no][col_no].red = triangles[i].clr.red;
                    frame_buffer[row_no][col_no].blue = triangles[i].clr.blue;
                    frame_buffer[row_no][col_no].green = triangles[i].clr.green;
                }
            }

        }
    }
//    for(int i = 0; i < screen_width; i++){
//        for(int j = 0; j < screen_height; j++){
//                //cout << z_buffer[i][j] << "\t";
//        }
//    }


    bitmap_image image(screen_width,screen_height);
    out4.open("z-buffer.txt");
    out4 << std::fixed;
    out4 << std::setprecision(6);
    for(int i = 0; i < screen_width; i++){
        for(int j = 0; j < screen_height; j++){
                if(z_buffer[i][j]<rear_z)
                {
                    out4 << z_buffer[i][j] << "\t" ;
                }
                image.set_pixel(j,i,frame_buffer[i][j].red,frame_buffer[i][j].green,frame_buffer[i][j].blue);
        }
        out4 << endl;
    }
    image.save_image("out.bmp");
    out4.close();

    for(int i = 0; i < screen_width; i++)
    {
        delete [] z_buffer[i];
        delete [] frame_buffer[i];
    }
    delete [] z_buffer;
    delete [] frame_buffer;
}



int main()
{

    string line;
    std::ifstream in("scene.txt");
    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to scene.txt

    //gluLookAt

    cin >> eye.x >> eye.y >> eye.z;
    cin >> look.x >> look.y >> look.z;
    cin >> up.x >> up.y >> up.z;

    //gluPerspective
    cin >> fov.y >> aspectRatio >> near >> far;

    view_transformation();
    projection_transformation();
    modeling_transformation();

    out.close();
    out2.close();
    out3.close();

    clip_scan_conversion();
}
