using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using System.Diagnostics;
using System.Threading;

namespace Histogram
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void btn1_Click(object sender, EventArgs e)
        {
            File.Delete("hist.csv");
            StreamWriter sw = new StreamWriter("hist.csv", true);
            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Filter = "BMP|*.bmp|JPG|*.jpg|PNG|*.png|ALL|*.*";
            if (ofd.ShowDialog() == DialogResult.OK)
            {
                if (ofd.FileName != "")
                {
                    Bitmap bmp = new Bitmap(ofd.FileName);
                    Bitmap bmp1 = new Bitmap(bmp.Width, bmp.Height);
                    pb1.Image = bmp;
                    int[] Ra = new int[256];
                    int[] Ga = new int[256];
                    int[] Ba = new int[256];
                    for (int i = 0; i < bmp.Width; i++)
                    {
                        for (int j = 0; j < bmp.Height; j++)
                        {
                            Color c = bmp.GetPixel(i, j);
                            Ra[c.R]++;
                            Ga[c.G]++;
                            Ba[c.B]++;
                            if (c.R > 150 && c.G > 150 && c.B < 100)
                                bmp1.SetPixel(i, j, Color.Blue);
                            else
                                bmp1.SetPixel(i, j, c);
                        }

                    }
                    pb2.Image = bmp1;
                    for (int x = 0; x < 256; x++)
                    {
                        sw.WriteLine(Ra[x].ToString() + "," + Ga[x].ToString() + "," + Ba[x].ToString());
                    }
                    sw.Close();
                }
            }
        }
        int count = 0;
        Bitmap bp1, bp2,bp3;
        private void t1_Tick(object sender, EventArgs e)
        {
            count = (count + 10) % 360;
            int x = 0, y = 0, minx = bp1.Width, maxx = 0, miny = bp1.Height,
            maxy = 0;
            int[,] xr = new int[bp1.Width, bp1.Height];
            int[,] yr = new int[bp1.Width, bp1.Height];
            Color[,] cc = new Color[bp1.Width, bp1.Height];
            for (int i = 0; i < bp1.Width; i++)
                for (int j = 0; j < bp1.Height; j++)
                {
                    x = Convert.ToInt32(i * Math.Cos(count * 3.14 / 180) + j * Math.Sin(count * 3.14 / 180));
                    y = Convert.ToInt32(j * Math.Cos(count * 3.14 / 180) - i * Math.Sin(count * 3.14 / 180));
                    xr[i, j] = x;
                    yr[i, j] = y;
                    cc[i, j] = bp1.GetPixel(i, j);
                    if (minx > x)
                        minx = x;
                    if (maxx < x)
                        maxx = x;
                    if (miny > y)
                        miny = y;
                    if (maxy < y)
                        maxy = y;
                }
            int fx = Math.Abs(minx) + 1;
            int fy = Math.Abs(miny) + 1;
            maxx += fx + 1;
            maxy += fy + 1;
            bp2 = new Bitmap(maxx, maxy);
            for (int i = 0; i < bp1.Width; i++)
                for (int j = 0; j < bp1.Height; j++)
                    bp2.SetPixel(xr[i, j] + fx, yr[i, j] + fy, cc[i, j]);
            pb2.Image = bp2;
        }

        private void btn_Rotation_Click(object sender, EventArgs e)
        {
            OpenFileDialog ofd = new OpenFileDialog();
            if (ofd.ShowDialog() == DialogResult.OK)
                if (ofd.FileName != "")
                {
                    pb1.Load(ofd.FileName);
                    bp1 = new Bitmap(pb1.Image);
                    bp2 = new Bitmap(bp1.Width * 2, bp1.Height * 2);
                    t1.Start();
                }
        }
        //Image Rotaion.
        private void btn_Stop_Click(object sender, EventArgs e)
        {
            t1.Stop();
        }

        private void btn_wave_Click(object sender, EventArgs e)
        {
            int iteration = Convert.ToInt32(tb1.Text);
            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Filter = "BMP|*.bmp|JPG|*.jpg|PNG|*.png|ALL|*.*";
            if (ofd.ShowDialog() == DialogResult.OK)
                if (ofd.FileName != "")
                {
                    pb1.Load(ofd.FileName);
                    bp1 = new Bitmap(pb1.Image);
                    int wid = bp1.Width;
                    int hgt = bp1.Height;
                    bp2 = new Bitmap(wid, hgt);
                    bp3 = new Bitmap(wid, hgt);
                    double[,] data = new double[wid, hgt];
                    for (int i = 0; i < wid; i++)
                    {
                        for (int j = 0; j < hgt; j++)
                        {
                            Color c = bp1.GetPixel(i, j);
                            data[i, j] = (c.R + c.G + c.B) / 3;
                        }

                    }
                    wavelet(data, iteration);
                    for (int i = 0; i < wid; i++)
                    {
                        for (int j = 0; j < hgt; j++)
                        {
                            int tmp= (int) data[i, j];
                            if (tmp > 255) tmp = 255;
                            if (tmp < 0) tmp = 0;
                            Color c = Color.FromArgb(tmp,tmp,tmp);
                            bp2.SetPixel(i, j, c);
                        }
                    }
                    bp2.Save("wave.jpg", System.Drawing.Imaging.ImageFormat.Jpeg);
                    pb2.Image = bp2;
                    Iwavelet(data, iteration);
                    for (int i = 0; i < wid; i++)
                    {
                        for (int j = 0; j < hgt; j++)
                        {
                            int tmp = (int)data[i, j];
                            if (tmp > 255) tmp = 255;
                            if (tmp < 0) tmp = 0;
                            Color c = Color.FromArgb(tmp, tmp, tmp);
                            bp3.SetPixel(i, j, c);
                        }
                    }
                    pb3.Image = bp3;
                }
        }
        //wavelet
        double s0 = 0.5, s1 = 0.5;
        double w0 = 1, w1 = -1;

        

        public void wavelet(double[] data)
        {
            double[] temp = new double[data.Length];
            int h = data.Length >> 1;
            for (int i = 0; i < h; i++)
            {
                int k = (i << 1);
                temp[i] =(data[k]);
                temp[i + h] = data[k] - data[k + 1];
            }
            for (int i = 0; i < data.Length; i++)
                data[i] = temp[i];
        }

        

        public void wavelet(double[,] data, int iterations)
        {
            int rows = data.GetLength(0);
            int cols = data.GetLength(1);
            double[,] data2 = new double[rows, cols]; 
            double[] row;
            double[] col;
            for (int k = 0; k < iterations; k++)
            {
                int lev = 1 << k;
                int levCols = cols / lev;
                int levRows = rows / lev;
                row = new double[levRows];
                for (int i = 0; i < levRows; i++)
                {
                    for (int j = 0; j < row.Length; j++)
                        row[j] = data[i, j];
                    wavelet(row);
                    for (int j = 0; j < row.Length; j++)
                        data[i, j] = row[j];
                }
                col = new double[levCols];
                for (int j = 0; j < levCols; j++)
                {
                    for (int i = 0; i < col.Length; i++)
                        col[i] = data[i, j];
                    wavelet(col);
                    for (int i = 0; i < col.Length; i++)
                        data[i, j] = col[i];
                }
                //data = data2;
            }
            
        }
        //inverse
        public void Iwavelet(double[] data)
        {
            double[] temp = new double[data.Length];
            int h = data.Length >> 1;
            for (int i = 0; i < h; i++)
            {
                int k = (i << 1);
                temp[k] = (data[i]);
                temp[k + 1] = (data[i]  + data[i + h] * w1);
            }
            for (int i = 0; i < data.Length; i++)
                data[i] = temp[i];
        }
        public void Iwavelet(double[,] data, int iterations)
        {
            int rows = data.GetLength(0);
            int cols = data.GetLength(1);
            double[] col;
            double[] row;
            for (int k = iterations - 1; k >= 0; k--)
            {
                int lev = 1 << k;
                int levCols = cols / lev;
                int levRows = rows / lev;
                col = new double[levRows];
                for (int j = 0; j < levCols; j++)
                {
                    for (int i = 0; i < col.Length; i++)
                        col[i] = data[i, j];
                    Iwavelet(col);
                    for (int i = 0; i < col.Length; i++)
                        data[i, j] = col[i];
                }
                row = new double[levCols];
                for (int i = 0; i < levRows; i++)
                {
                    for (int j = 0; j < row.Length; j++)
                        row[j] = data[i, j];
                    Iwavelet(row);
                    for (int j = 0; j < row.Length; j++)
                        data[i, j] = row[j];
                }
            }
        }
        //////////////////////////////////////////////////////////////
        //////////////////////Sobel Filter////////////////////////////
        //////////////////////////////////////////////////////////////
        
        private void btn_sobel_Click(object sender, EventArgs e)
        {
            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Filter = "BMP|*.bmp|JPG|*.jpg|PNG|*.png|ALL|*.*";
            if (ofd.ShowDialog() == DialogResult.OK)
                if (ofd.FileName != "")
                {
                    pb1.Load(ofd.FileName);
                    bp1 = new Bitmap(pb1.Image);
                    int wid = bp1.Width;
                    int hgt = bp1.Height;
                    bp2 = new Bitmap(wid, hgt);
                    double[,] data = new double[wid, hgt]; // gray level
                    for (int i = 0; i < wid; i++)
                    {
                        for (int j = 0; j < hgt; j++)
                        {
                            Color c = bp1.GetPixel(i, j);
                            data[i, j] = (c.R + c.G + c.B) / 3;
                        }

                    }
                    double gx;
                    double gy;
                    double val;
                    double[,] dest = new double[wid, hgt];
                    for (int x = 1; x < wid - 1; x++)
                    {
                        for (int y = 1; y < hgt - 1; y++)
                        {
                            //                     sobel-y --> [1   2   1]
                            //                                 [0   0   0]
                            //                                 [-1 -2  -1]

                            //                    sobel-x -->  [1 0 -1]
                            //                                 [2 0 -2]
                            //                                 [1 0 -1] 

                            gx = data[x - 1, y - 1] + 2 * data[x, y - 1] + data[x + 1, y - 1] - data[x - 1, y + 1] - 2 * data[x, y + 1] - data[x + 1, y + 1];
                            gy = data[x - 1, y - 1] + 2 * data[x - 1, y] + data[x - 1, y + 1] - data[x + 1, y - 1] - 2 * data[x + 1, y] - data[x + 1, y + 1];
                            val = Math.Abs(gx) + Math.Abs(gy);
                            dest[x, y] = val;
                        }
                    }
                    data = dest;
                    for (int i = 0; i < wid; i++)
                    {
                        for (int j = 0; j < hgt; j++)
                        {
                            int tmp = (int)data[i, j];
                            if (tmp > 255) tmp = 255;
                            if (tmp < 0) tmp = 0;
                            Color c = Color.FromArgb(tmp, tmp, tmp);
                            bp2.SetPixel(i, j, c);
                        }
                    }
                    pb2.Image = bp2;
                }
        }

        

        ///////////////////////////////////////////////////
        /// ///////////Seed Filling//////////////////////////////
        /// ///////////////////////////////////////////
        public void seedfill(Bitmap bmp,Color c,Point point)
        {
            //Remove the color of the seed point (a point within the boundary of the graphic that needs to be filled)
            bool[,] flg=new bool[bmp.Width, bmp.Height];
           List<Point> points = new List<Point>();
            bmp.SetPixel(point.X, point.Y, Color.Red);
            points.Add(point);
            while (points.Count > 0)
            {
                if ( point.X + 1 < bmp.Width && bmp.GetPixel(point.X + 1, point.Y) == c && flg[point.X + 1, point.Y] == false && !points.Exists(i=>i == (new Point(point.X+1, point.Y))))
                { points.Add(new Point(point.X + 1, point.Y)); }
                if (point.Y +1 < bmp.Height && bmp.GetPixel(point.X , point.Y+1) == c && flg[point.X, point.Y+1] == false && !points.Exists(i => i == (new Point(point.X , point.Y + 1))))
                { points.Add(new Point(point.X , point.Y + 1)); }
                if (point.X - 1 > 0 && bmp.GetPixel(point.X - 1, point.Y) == c && flg[point.X - 1, point.Y] == false && !points.Exists(i => i == (new Point(point.X - 1, point.Y))))
                { points.Add(new Point(point.X - 1, point.Y)); }
                if (point.Y - 1  > 0 && bmp.GetPixel(point.X , point.Y-1) == c && flg[point.X , point.Y-1] == false && !points.Exists(i => i == (new Point(point.X , point.Y - 1))))
                { points.Add(new Point(point.X , point.Y-1)); }
                point = points.Last();
                points.Remove(point);
                flg[point.X, point.Y] = true;
                bmp.SetPixel(point.X, point.Y, Color.Red);
                points = points.Distinct().ToList();
                pb2.Refresh();
            }
        }

        

        bool[,] flg;

        
        public void seedfill2(Bitmap bmp, Color c, Point point)
        {
            //Remove the color of the seed point (a point within the boundary of the graphic that needs to be filled)
            bmp.SetPixel(point.X, point.Y, Color.Red);
            flg[point.X, point.Y] = true;
            //pb2.Refresh();
            if (point.X + 1 < bmp.Width && bmp.GetPixel(point.X + 1, point.Y) == c && flg[point.X + 1, point.Y] == false)
            { seedfill2(bmp, c, new Point(point.X + 1, point.Y)); }
            if (point.Y + 1 < bmp.Height && bmp.GetPixel(point.X, point.Y + 1) == c && flg[point.X, point.Y + 1] == false )
            { seedfill2(bmp, c, new Point(point.X, point.Y + 1)); }
            if (point.X - 1 > 0 && bmp.GetPixel(point.X - 1, point.Y) == c && flg[point.X - 1, point.Y] == false )
            { seedfill2(bmp, c, new Point(point.X - 1, point.Y)); }
            if (point.Y - 1 > 0 && bmp.GetPixel(point.X, point.Y - 1) == c && flg[point.X, point.Y - 1] == false )
            { seedfill2(bmp, c, new Point(point.X, point.Y - 1)); }

        }
        private void btn_seed_Click(object sender, EventArgs e)
        {
            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Filter = "BMP|*.bmp|JPG|*.jpg|PNG|*.png|ALL|*.*";
            
            if (ofd.ShowDialog() == DialogResult.OK)
                if (ofd.FileName != "")
                {
                    pb1.Load(ofd.FileName);
                    bp1 = new Bitmap(ofd.FileName);
                    int wid = bp1.Width;
                    int hgt = bp1.Height;
                    bp2 = new Bitmap(wid, hgt);
                    flg=new bool[wid,hgt];
                    pb2.Image = bp1;
                    Color c= bp1.GetPixel(wid/2, hgt/2);
                    Thread th1 = new Thread( ()=> { seedfill2(bp1, c, new Point(wid / 2, hgt / 2)); });
                    Thread th2 = new Thread(() => { seedfill2(bp1, c, new Point(wid / 4, hgt / 4)); });
                    th1.Start();
                    th2.Start();
                    //seedfill(bp1,c,new Point(wid/2, hgt/2));
                    
                }
        }
        // /////////////////////////////////////////////////////////////////////
        // /////////////////////////shuffle////////////////////////////////////
        // ///////////////////////////////////////////////////////////////////
        private void btn_shuffle_Click(object sender, EventArgs e)
        {
            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Filter = "ALL|*.*|BMP |*.bmp|JPG|*.jpg|PNG|*.png";
            if (ofd.ShowDialog() == DialogResult.OK)
            {
                if (ofd.FileName != "")
                {
                    Bitmap bmp = new Bitmap(ofd.FileName);
                    Bitmap bmp1 = new Bitmap(bmp.Width, bmp.Height);
                    pb1.Image = bmp;
                    int wid = bmp.Width; int hgt = bmp.Height;
                    for (int i = 0; i < (wid / 2) - 1; i++)
                    {
                        for (int j = 0; j < (hgt / 2) - 1; j++)
                        {
                            Color c = bmp.GetPixel(i, j);
                            bmp1.SetPixel(i + (wid / 2), j + (hgt / 2), c);

                            Color c1 = bmp.GetPixel(i + (wid / 2), j + (hgt / 2));
                            bmp1.SetPixel(i, j, c1);
                        }
                    }
                    for (int i = (wid / 2); i < wid - 1; i++)
                    {
                        for (int j = 0; j < (hgt / 2) - 1; j++)
                        {
                            Color c = bmp.GetPixel(i, j);
                            bmp1.SetPixel(i - (wid / 2), j + (hgt / 2), c);

                            Color c1 = bmp.GetPixel(i - (wid / 2), j + (hgt / 2));
                            bmp1.SetPixel(i, j, c1);
                        }
                    }
                    pb2.Image = bmp1;
                }
            }
            
        }
        // /////////////////////////////////////////////////////////////////////
        // /////////////////////////sqrRGB////////////////////////////////////
        // ///////////////////////////////////////////////////////////////////
        private void btn_sqrRGB_Click(object sender, EventArgs e)
        {
            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Filter = "ALL|*.*|BMP |*.bmp|JPG|*.jpg|PNG|*.png";
            if (ofd.ShowDialog() == DialogResult.OK)
            {
                if (ofd.FileName != "")
                {
                    Bitmap bmp = new Bitmap(ofd.FileName);
                    Bitmap bmp1 = new Bitmap(bmp.Width, bmp.Height);
                    pb1.Image = bmp;
                    int wid = bmp.Width; int hgt = bmp.Height;
                    for (int i = 0; i < wid - 1; i+=2)
                    {
                        for (int j = 0; j < hgt - 1; j+=2)
                        {
                            Color c = bmp.GetPixel(i, j);
                            bmp1.SetPixel(i /2, j /2, c);
                            bmp1.SetPixel((i/2)+(wid/2), j/2, Color.FromArgb(c.R,0,0));
                            bmp1.SetPixel((i / 2) , (j/2)+(hgt/2), Color.FromArgb(0, c.G, 0));
                            bmp1.SetPixel((i / 2) + (wid / 2), (j / 2) + (hgt / 2), Color.FromArgb(0, 0, c.B));
                        }
                    }
                    pb2.Image = bmp1;
                }
            }
        }
        // /////////////////////////////////////////////////////////////////////
        // /////////////////////////halfRGB////////////////////////////////////
        // ///////////////////////////////////////////////////////////////////
        private void btn_RGB_Click(object sender, EventArgs e)
        {
            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Filter = "ALL|*.*|BMP |*.bmp|JPG|*.jpg|PNG|*.png";
            if (ofd.ShowDialog() == DialogResult.OK)
            {
                if (ofd.FileName != "")
                {
                    Bitmap bmp = new Bitmap(ofd.FileName);
                    Bitmap bmp1 = new Bitmap(bmp.Width, bmp.Height);
                    Bitmap bmp2 = new Bitmap(bmp.Width, bmp.Height);
                    Bitmap bmp3 = new Bitmap(bmp.Width, bmp.Height);
                    pb1.Image = bmp;
                    int wid = bmp.Width; int hgt = bmp.Height;
                    for (int i = 0; i < wid ; i++)
                    {
                        for (int j = 0; j < hgt ; j++)
                        {
                            Color c = bmp.GetPixel(i, j);
                            Color r = Color.FromArgb(c.R, 0, 0);
                            Color g = Color.FromArgb(0, c.G, 0);
                            Color b = Color.FromArgb(0, 0, c.B);
                            if (i > wid / 2)
                            {
                                bmp1.SetPixel(i, j, c);
                                bmp2.SetPixel(i, j, c);
                                bmp3.SetPixel(i, j, c);
                            }
                            else
                            {
                                bmp1.SetPixel(i, j, r);
                                bmp2.SetPixel(i, j, g);
                                bmp3.SetPixel(i, j, b);
                            }
                        }
                    }
                    pb2.Image = bmp1;
                    pb3.Image = bmp2;
                    pb4.Image = bmp3;
                }
            }
        }

    }
}
