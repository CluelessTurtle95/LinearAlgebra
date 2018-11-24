// Program to calculate the Determinant,Inverse,Eigenvalues(hopefully) of any square matrix
// Multiply matrices together, by scalars , add them etc.

#include "pch.h"
#include<iostream>
#include<stdlib.h>
#include<cmath>
#include<conio.h>
#include<math.h>
#include<fstream>
#include<iomanip>

using namespace std;
#define SIZEc 3 // default conversion

int SIZE = SIZEc;
/*
void createspace();
void openspace();
void createtrans();
void composetrans();
void edittrans();
void vectoralg();
*/
float approx(float a)
{
	if (abs(a) < 0.000001)
		a = 0;

	return a;
}

class matrix
{
	// :private
	int size;
	float* *m; //= new float*[size]; // create a dynamically allocated array of pointers 
								   //each of which pointing to another dynamically allocated array = 2D ARRAY 
	float trace;
	float det;
	int id;

	void createid()
	{
		fstream file;
		file.open("maDATA.dat", ios::in | ios::binary);

		if (file.is_open())
		{
			file.read((char *)this, sizeof(matrix));
			id++;
		}
		else
		{
			id = 1;
		}
	}

	float cofactor(int i, int j)
	{
		if (size == 1)
		{
			return m[0][0];
		}

		float ** a = new float *[size - 1];

		for (int i = 0; i < size - 1; i++)
		{
			a[i] = new float[size - 1];
		}

		int p, q;

		for (int l = 0; l < size; l++)
			for (int n = 0; n < size; n++)
			{
				p = l;
				q = n;
				if (l != i && n != j)
				{
					if (l > i)
					{
						p--;
					}
					if (n > j)
					{
						q--;
					}

					a[p][q] = m[l][n];
				}
			}

		matrix t(a, size - 1);

		return pow(-1, (i + j)) * t.getdet();
	}

public:

	matrix(float **a, int s)
	{
		size = s;
		m = new float*[size];

		for (int i = 0; i < size; i++)
		{
			m[i] = new float[size];
		}

		m = a;

		ftrace();
		getdet();
		createid();
		save();
	}

	int getid()
	{
		return id;
	}

	void save()
	{
		fstream infile;
		infile.open("maDATA.dat", ios::in | ios::binary);

		matrix a(m , size);
		int pos = 0;
		while (!infile.eof())
		{
			infile.read((char *)&a, sizeof(matrix));
			if (a.getid() == id)
			{
				pos = infile.tellg();
			}
		}
		if (pos != 0)
			pos -= sizeof(matrix);

		infile.open("maDATA.dat", ios::out | ios::binary);
		infile.seekp(pos);
		infile.write((char *)this, sizeof(matrix));
	}

	/*
	~matrix()
	{
		for(int i=0;i < size;i++)
		delete m[i];

		delete m;
	}
	*/ // destructor

	void clean()
	{
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			{
				m[i][j] = approx(m[i][j]);
			}
	}

	void ftrace()
	{
		trace = 0;

		for (int i = 0; i < size; i++)
		{
			trace += m[i][i];
		}
	}

	float gettrace()
	{
		return trace;
	}

	int getsize()
	{
		return size;
	}

	float element(int i, int j)
	{
		return m[i][j];
	}

	float getdet()
	{
		det = 0;

		if (size == 1)
		{
			det = m[0][0];
			return det;
		}

		for (int i = 0; i < size; i++)
		{
			det += m[0][i] * cofactor(0, i);
		}

		clean();

		return det;
	}

	float ** cofactormatrix()
	{
		float ** cm = new float*[size];

		for (int i = 0; i < size; i++)
		{
			cm[i] = new float[size];
		}

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				cm[i][j] = cofactor(i, j);
			}
		}

		return cm;
	}

	void transpose()
	{
		float temp;

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			{
				if (i > j)
				{
					temp = m[i][j];
					m[i][j] = m[j][i];
					m[j][i] = temp;
				}
			}
		save();
	}

	void inverse()
	{
		if (det == 0)
		{
			cout << " NOT INVERTIBLE! " << endl;
		}

		m = cofactormatrix();
		transpose();

		if (det)
			for (int i = 0; i < size; i++)
				for (int j = 0; j < size; j++)
				{
					m[i][j] = m[i][j] / det;
				}
		clean();
		ftrace();
		save();
	}

	static float** convert(float a[SIZEc][SIZEc], int s)
	{
		float **x = new float*[s];

		for (int i = 0; i < s; i++)
		{
			x[i] = new float[s];
		}

		for (int j = 0; j < s; j++)
			for (int i = 0; i < s; i++)
			{
				x[i][j] = a[i][j];
			}

		return x;
	}

	void print()
	{
		clean();
		for (int i = 0; i < size; i++)
		{
			cout << endl;
			for (int j = 0; j < size; j++)
				cout << m[i][j] << " ";
		}
	}

	void modify(int i, int j, float value)
	{
		m[i][j] = value;
		ftrace();
		save();
	}

	void scalar(float value)
	{
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			{
				m[i][j] = value * m[i][j];
			}
		ftrace();
		save();
	}

	static matrix multi(matrix a, matrix b)
	{
		float t = 0;

		matrix c = matrix::null();

		for (int i = 0; i < a.getsize(); i++)
		{
			t = 0;

			for (int j = 0; j < a.getsize(); j++)
			{
				for (int k = 0; k < a.getsize(); k++)
				{
					t += a.element(k, i) * b.element(j, k);
				}
				c.modify(i, j, t);
				t = 0;
			}
		}
		return c;
	}

	static matrix add(matrix a, matrix b)
	{
		if (a.getsize() == b.getsize())
		{
			matrix c(a.cofactormatrix(), a.getsize()); // just random initializasion
			for (int i = 0; i < a.getsize(); i++)
				for (int j = 0; j < a.getsize(); j++)
					c.modify(i, j, a.element(i, j) + b.element(i, j));
			return c;
		}
		else
			return matrix::null();
	}

	static matrix null()
	{
		float **a = new float*[SIZE];
		for (int i = 0; i < SIZE; i++)
			a[i] = new float[SIZE];

		matrix temp(a, SIZE);
		return temp;
	}
};

class vector
{
	int size;
	float * v;
	int flag;
	int id;
	float mag;

	void createid()
	{
		fstream file;
		file.open("vecDATA.dat", ios::in | ios::binary);
		
		if (file.is_open())
		{
			file.read((char *)this, sizeof(vector));
			id++;
		}
		else
		{
			id = 1;
		}
	}

public:
	int getid() 
	{
		return id;
	}

	vector(float *x, int s)
	{
		size = s;
		v = new float[size];
		v = x;
		flag = 1;
		createid();
		cout << "AAA";
		save();
	}

	vector()
	{
		flag = 0;
	}

	float element(int i)
	{
		return v[i];
	}

	int getsize()
	{
		return size;
	}

	void save() 
	{
		fstream infile;
		infile.open("vecDATA.dat" , ios::in|ios::binary);
		int pos = 0;
		if (infile.is_open())
		{
			vector a;
			
			while (!infile.eof())
			{
				infile.read((char *)&a, sizeof(vector));
				if (a.getid() == id)
				{
					pos = infile.tellg();
				}
			}
			if (pos != 0)
				pos -= sizeof(vector);
		}

		infile.open("vecDATA.dat", ios::out | ios::binary);
		infile.seekp(pos);
		infile.write((char *)this, sizeof(vector));
	}
	
	void modify(int i, float value)
	{
		v[i] = value;
		save();
	}

	static vector crossproduct(vector a, vector b)
	{
		vector c = a;

		int g = c.getsize();

		float ** d = new float *[g];

		for (int i = 0; i < g; i++)
		{
			d[i] = new float[g];
		}


		for (int k = 0; k < g; k++)
		{
			for (int i = 0; i < g; i++)
			{
				for (int j = 0; j < g; j++)
				{
					if (i == 1)
					{
						d[i][j] = a.element(j);
					}
					if (i == 2)
					{
						d[i][j] = b.element(j);
					}
				}
			}
			for (int i = 0; i < g; i++)
				if (i == k)
				{
					d[0][i] = 1;
				}
				else
				{
					d[0][i] = 0;
				}

			matrix e(d, g);
			c.modify(k, e.getdet());
		}
		return c;
	}

	static float dotproduct(vector a, vector b)
	{
		float c = 0;

		if (a.getsize() != b.getsize())
		{
			return -1;
		}

		for (int i = 0; i < a.getsize(); i++)
		{
			c += a.element(i)*b.element(i);
		}

		return c;
	}

	static vector null(int s)
	{

		float * a = new float[s];
		vector c = vector(a, s);
		return c;
	}

	void scalar(float k)
	{
		for (int i = 0; i < size; i++)
		{
			v[i] = k * v[i];
		}
		save();
	}

	static vector add(vector a, vector b)
	{
		vector c = a;
		float t = 0;
		for (int i = 0; i < a.getsize(); i++)
		{
			t = a.element(i) + b.element(i);
			c.modify(i, t);
		}
		return c;
	}

	void print()
	{
		for (int i = 0; i < size; i++)
		{
			cout << v[i] << endl;
		}
	}

	void magnitude()
	{
		mag = 0;
		for (int i = 0; i < size; i++)
		{
			mag += pow(v[i], 2);
		}
		mag = pow(mag, 0.5);
	}

	float getmag()
	{
		return mag;
	}

	static bool compare(vector a, vector b)
{
	bool result = true;
	if (a.getsize() == b.getsize())
	{
		for (int i = 0; i < a.getsize(); i++)
		{
			if (a.element(i) != b.element(i))
				result = false;
		}
	}
	else
	{
		result = false;
	}

	return result;
}
};

class space
{
	vector * s;
	int elements;
	int dimensions;
	int count;

	int id;

	void createid()
	{
		fstream file;
		file.open("spaDATA.dat", ios::in | ios::binary);

		if (file.is_open())
		{
			file.read((char *)this, sizeof(space));
			id++;
		}
		else
		{
			id = 1;
		}
	}
public:

	space(int x, int size)
	{
		elements = x;
		count = 0;

		dimensions = size;
		s = new vector[elements];

		float *a = new float[dimensions];

		for (int i = 0; i < dimensions; i++)
			a[i] = 0;

		for (int i = 0; i < elements; i++)
			s[i] = vector(a, dimensions);

		createid();
	}

	int getid()
	{
		return id;
	}

	void save()
	{
		update();
		fstream infile;
		infile.open("spaDATA.dat", ios::in | ios::binary);

		space a = *this;
		int pos = 0;
		while (!infile.eof())
		{
			infile.read((char *)&a, sizeof(space));
			if (a.getid() == id)
			{
				pos = infile.tellg();
			}
		}
		if (pos != 0)
			pos -= sizeof(space);

		infile.open("spaDATA.dat", ios::out | ios::binary);
		infile.seekp(pos);
		infile.write((char *)this, sizeof(space));
	}

	void add(vector a)
	{
		if (a.getsize() != dimensions)
		{
			return;
		}

		if (count < elements)
		{
			s[count] = a;
			count++;
		}
		else
		{
			return;
		}
		save();
	}

	void del(int m)
	{
		if (count > m)
		{
			for (int i = m; i < count - 1; i++)
			{
				s[i] = s[i + 1];
			}
		}
		else
		{
			count++;
		}
		count--;
		save();
	}

	void del(vector a)
	{
		for (int i = 0; i < count; i++)
			if (vector::compare(s[i], a))
			{
				del(i);
			}
	}

	vector at(int i)
	{
		vector temp = vector::null(SIZE);

		if (i < count)
			return s[i];
		else
			return temp;
	}

	void apply(matrix a)
	{
		for (int m = 0; m < count; m++)
		{
			vector c = s[m].null(a.getsize());

			if (a.getsize() != s[m].getsize())
			{
				return;
			}

			float  t;

			for (int i = 0; i < s[m].getsize(); i++)
			{
				t = 0;
				for (int k = 0; k < a.getsize(); k++)
				{
					t += a.element(i, k)*s[m].element(k);
				}

				c.modify(i, t);
			}
			save();
		}

	}

	int getelements()
	{
		return elements;
	}

	int getdimensions()
	{
		return dimensions;
	}

	void print()
	{
		cout << setw(4) << "ID " << " | " << setw(9) << "Vector";
		for (int i = 0; i < count; i++)
		{
			cout << setw(4) << s[i].getid() << " | ";
			for (int i = 0; i < SIZE; i++)
				cout << setw(3) << s[i].element(i);
		}
	}

	void update()
	{
		fstream infile;
		infile.open("vecDATA.dat" , ios::binary|ios::in);
		if (!infile.is_open())
		{
			return;
		}

		vector h;
		for (int i = 0; i < count; i++)
		{
			infile.seekg((s[i].getid() - 1) * sizeof(vector), ios::beg);
			infile.read((char*)&h, sizeof(h));
			s[i] = h;
		}
		infile.close();
	}
};

/*
class graph
{
	int x, y, z;
	space s;
public:
	graph(space sp)
	{
		s = sp;
	}
};
*/

static class menu
{
public:
	static void createspace()
	{
		fstream file;
		file.open("vecDATA.dat" , ios::binary | ios::in);
		
		if (!file.is_open())
		{
			cout << "NO VECTORS AVAILABLE!" << endl;
			return;
		}

		vector h;
		int n ;
		
		cout << "Number of Vectors : ";
		cin >> n;

		int *a = new int[n];

		cout << "Chose vectors to add:" << endl;

		cout << setw(4) << "ID " << " | " << setw(9) << "Vector";

		while (!file.eof())
		{
			file.read((char*)&h, sizeof(vector));
			cout << setw(4)<< h.getid() << " | ";
			for (int i = 0; i < SIZE; i++)
				cout << setw(3) << h.element(i);
		}
		int vecnum = h.getid();

		cout << endl << "Chose Vectors : ";
		for (int i = 0; i < n; i++)
		{
			cin >> a[i];
		}
		space  s( n  , SIZE);

		file.seekg(0, ios::beg);

		int i = 0;
		while (!file.eof())
		{
			file.read((char*)&h, sizeof(vector));
			if (a[i] == h.getid())
				s.add(h);
			i++;
		}
		s.save();

		cout << "SPACE ID IS : " << s.getid();
	}

	static void openspace()
	{
		cout << "Enter Space ID : ";
		int f, e;
		cin >> f;

		fstream infile;
		infile.open("spaDATA.dat", ios::in | ios::binary);

		if (!infile.is_open())
		{
			cout << "NO SPACES FOUND!" << endl;
			infile.close();
			return;
		}

		space a = space(2, 2);

		while (!infile.eof())
		{
			infile.read((char*)&a, sizeof(space));
			if (a.getid() == f)
			{
				break;
			}
		}
		cout << "Space " << a.getid() << " Opened !" << endl;

		cout << "1. Apply Transformation \n";
		cout << "2. Edit Space\n";
		cout << "3. Back\n";

		cout << "Choice Number : ";
		cin >> f;

		switch (f)
		{
		case 1:
			{
				infile.open("maDATA.dat", ios::in | ios::binary);
					
				if (!infile.is_open())
				{
					cout << "NO Available Transformation!";
					break;
				}

				cout << " Enter Transformation ID : ";
				cin >> e;
				float **b = new float*[SIZE];
				matrix m = matrix(b, e);

				while (!infile.eof())
				{
					infile.read((char*)&m, sizeof(matrix));

					if (m.getid() == e)
					{
						a.apply(m);
						cout << "APPLIED! " << endl;
						break;
					}
				}
				break;
			}
		case 2:
			{
				cout << endl << " Current Space : " << endl;
				a.print();

				cout << "1. Add Vector \n"
					<< "2. Remove Vector \n"
					<< "3. Back\n";
	
				int b;
	
				cout << "Choice Number : ";
				cin >> b;

				switch (b)
				{
				case 1:
					{
						int c;
						cout << "Enter Vector ID :";
						cin >> c;
		
						fstream infile("vecDATA.dat", ios::in | ios::binary);
						if (!infile.is_open())
						{
							cout << "NO VECTORS FOUND!";
						}

						vector h;

						while (!infile.eof())
						{
							infile.read((char*)&h, sizeof(vector));
							if (c == h.getid())
							{
								a.add(h);
								cout << "Added!";
							}
						}
						break;
					}

				case 2:
					{
						int c;
						cout << "Enter Vector ID :";
						cin >> c;

						fstream infile("vecDATA.dat", ios::in | ios::binary);
						if (!infile.is_open())
						{
							cout << "NO VECTORS FOUND!";
						}

						vector h;

						while (!infile.eof())
						{
							infile.read((char*)&h, sizeof(vector));
							if (c == h.getid())
							{
								a.del(h);
								cout << "Removed!";
							}
						}
						break;
					}

				case 3:
					break;
				}

			break;
			}

		case 3:
			break;
		}
	}
	
	static void createtrans()
	{

	}

	static void composetrans()
	{
	}

	static void edittrans()
	{

	}

	static void vectoralg()
	{
		cout << "\nVECTOR MENU \n\n";
		cout << "1. Add New Vector\n";
		cout << "2. Edit a Vector\n";
		cout << "3. Dot Product\n";
		cout << "4. Cross Product \n";
		cout << "5. Comparisons \n";
		cout << "6. Magnitude\n";
		
		int a,b;
		cout << "Choice Number : ";
		cin >> a;

		switch (a)
		{
			case 1:
			{
				cout << "SIZE : ";
				cin >> b;
			
				float *x = new float[b];
		
				cout << "Enter Elements : ";
				for (int i = 0; i < b; i++)
					cin >> x[i];
				
				cout << "dane q";
				vector v(x, b);
				cout << "done";
				break;
			}
				
			case 2:
			{
				float v;
				cout << "ID : ";
				cin >> b;

				fstream file("vecDATA.dat", ios::in | ios::binary);
	
				if (!file.is_open())
				{
					file.close();
					break;
				}	
					
				vector h;
				
				while (!file.eof())
				{
					file.read((char *)&h, sizeof(vector));
					if (h.getid() == b)
					{
						cout << "Vector Found!\n";
						h.print();
						cout << "\nElement No.(0-" << h.getsize() - 1 << ") : ";
						cin >> b;
						cout << "Value : ";
						cin >> v;
						h.modify(b, v);
					}

					break;
				}
			}
	
			case 3:
				break;
			case 4:
				break;
			case 5:
				break;
			case 6:
				break;
		}
	}
};

void mainmenu()
{
	system("clear");

	cout << "**********LINEAR-ALGEBRA**********" << endl;
	cout << "MAIN MENU" << endl;

	cout << "1. Create New Space" << endl;
	cout << "2. Open Existing Space" << endl;
	cout << "3. Create New Transformation" << endl;
	cout << "4. Compose Transformations" << endl;
	cout << "5. Edit transformations" << endl;
	cout << "6. Vector algebra" << endl;

	int choice;
	cout << "Choice Number : ";
	cin >> choice;

	switch (choice)
	{
	case 1:
		menu::createspace();
		break;
	case 2:
		menu::openspace();
		break;
	case 3:
		menu::createtrans();
		break;
	case 4:
		menu::composetrans();
		break;
	case 5:
		menu::edittrans();
		break;
	case 6:
		menu::vectoralg();
		break;
	case 7:
		exit(0);
		break;

	default:
		cout << "INVALID CHOICE!";
	}
}

int main()
{
	// const int SIZE = stoi(argv[1]);
	cin >> SIZE;
	mainmenu();
	_getch();
	// space s; s.create, 
	return 0;
}
