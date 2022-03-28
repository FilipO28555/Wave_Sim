#define OLC_PGE_APPLICATION
#define _USE_MATH_DEFINES
#include<Filip/olcPixelGameEngine.h>
#include<Filip/LittleHelp.h>



const int SIZEX = 1920;
const int SIZEY = 1080;

const int NX = 1920/5 ;
const int NY = NX*SIZEY/SIZEX ;

bool QuadratureVelocity = false;
bool QuadraturePosition = false;

double dt = 0.5;  // timestep
double h = 1;	  // spacestep

using olc::vd2d;
using std::cout;

struct cell {
	double dx[2],v,c;		// deviation/position (for fast euler we need two), velocity, speed of light in cell
	double* cdx = &dx[0];  //currently used dx
	cell* N[4];			  //	 Neughbours
	double grad[2];		 // = (gradx, grady). left - right, up - down.
	bool active = true;


	//Quadratic extrapolation:
	double a1, a2;
	double v1, v2;
	cell() {
	}
	cell(double DX, double C) {
		a1 = 0;
		a2 = 0;
		v1 = 0;
		v2 = 0;
		dx[0] = DX;
		dx[1] = DX;
		c = C;
		v = 0;
		for (int i = 0; i < 4; i++)
			N[i] = 0;
	}
	cell(double DX, double C, double V) {
		a1 = 0;
		a2 = 0;
		v1 = 0;
		v2 = 0;
		dx[0] = DX;
		dx[1] = DX;
		c = C;
		v = V;
		for (int i = 0; i < 4; i++)
			N[i] = 0;
	}
	void clearQuad() {
		a1 = 0;
		a2 = 0;
		v1 = 0;
		v2 = 0;
	}
	//dt dependent on elapsed frame time. Using quadratic extrapolation assuming dt does not change.
	void updateVQ() {
		double DX = 0;
		for (int i = 0; i < 4; i++) {
			//std::cout << N[i] << "\n";
			if (N[i]) {
				double temp = (N[i]->dx[0] - dx[0]);
				DX += temp;
				grad[i / 2] = (1 - 2 * (i % 2)) * temp;
			}
			else {
				double temp = v*dt * h / c; //Dumpening on edges
				DX -= temp;
				grad[i / 2] = (1 - 2 * (i % 2)) * temp;
			}
		}
		double a3 = (DX)/(h*h)*(c * c);
		v += (3*a3-3*a2+a1) * dt;

		a1 = a2;
		a2 = a3;
	}
	//Quad also
	void updateDxQ() {
		dx[0] += (3 * v - 3 * v2 + v1) * dt;
		v1 = v2;
		v2 = v;

	}
	void updateVEuler() {
		double DX = 0;
		for (int i = 0; i < 4; i++)
			if (N[i] && N[i]->active)
				DX += (N[i]->dx[0] - dx[0]);
			else {
				DX -= v * h; //Dumpening on edges
			}
		double a = c * c * (DX) / (h * h);
		v += a * dt;

	}
	void updateDxEuler() {
		dx[0] += v*dt;
	}

	void updateFastEuler(int par) {
		double DX = 0;
		for (int i = 0; i < 4; i++)
			if (N[i] && N[i]->active)
				DX += (N[i]->dx[par]);
			else {
				DX -= (h / dt - 1) * dx[par] - dx[1 - par] * h / dt; //Dumpening on edges
			}
		double y2 = c * c * dt * dt / h / h;
		dx[1 - par] = 2 * dx[par] - dx[1 - par] + (dt * dt * c * c / h / h) * (DX - 4 * dx[par]);
		//cdx = &dx[1 - par];
	}


};

struct Tabcell {
	double dx[2];
	double c; //prêdkoœæ œwiat³a
	Tabcell* N[4];   //s¹siedzi
	bool active = true;


	//Quadratic extrapolation:
	double a1, a2;
	double v1, v2;
	Tabcell() {
	}
	Tabcell(double DX, double C) {
		a1 = 0;
		a2 = 0;
		v1 = 0;
		v2 = 0;
		dx[0]=dx[1] = DX;
		c = C;
		
		for (int i = 0; i < 4; i++)
			N[i] = 0;
	}
	Tabcell(double DX, double C, double V) {
		a1 = 0;
		a2 = 0;
		v1 = 0;
		v2 = 0;
		dx[0] = dx[1] = DX;
		c = C;
		for (int i = 0; i < 4; i++)
			N[i] = 0;
	}
	//dt dependent on elapsed frame time. Using quadratic extrapolation assuming dt does not change.
	void update(int par) {
		double DX = 0;
		for (int i = 0; i < 4; i++)
			if (N[i] && N[i]->active)
				DX += (N[i]->dx[par] - dx[par]);
			else {
				DX -= (dx[par]-dx[1-par])*dt*h * c; //Dumpening on edges
			}
		double a3 = (DX) / (h * h) * (c * c);
		dx[1 - par] -= ((3 * a3 - 3 * a2 + a1)*dt) * dt + dx[1-par] ;

		a1 = a2;
		a2 = a3;
	}
};

std::array<std::array<cell, NY>, NX> Field;

template <typename T>
void Join(std::array<std::array<T, NY>, NX> &X) {
	for (int i = 0; i < NX; i++)
		for (int j = 0; j < NY; j++) {
			if (i + 1 < NX)
				X[i][j].N[0] = &X[i + 1][j]; //N[0] -> prawy
			else
				X[i][j].N[0] = 0;
			if (i - 1 >= 0)
				X[i][j].N[1] = &X[i - 1][j]; //N[1] -> lewy
			else
				X[i][j].N[1] = 0;
			if (j + 1 < NY)
				X[i][j].N[2] = &X[i][j + 1]; //N[2] -> dolny
			else
				X[i][j].N[2] = 0;
			if (j - 1 >= 0)
				X[i][j].N[3] = &X[i][j - 1]; //N[3] -> górny
			else
				X[i][j].N[3] = 0;
		}
}

std::chrono::duration<double> elapsedTime;
void Update(int parity) {
	// START TIMING
	auto tp1 = std::chrono::high_resolution_clock::now();
	if (QuadratureVelocity || QuadraturePosition) {
		if (QuadratureVelocity) {
			for (int i = 0; i < NX; i++)
				for (int j = 0; j < NY; j++) {
					Field[i][j].updateVQ();
				}
		}
		else {
			for (int i = 0; i < NX; i++)
				for (int j = 0; j < NY; j++) {
					Field[i][j].updateVEuler();
				}
		}
		if (QuadraturePosition) {
			for (int i = 0; i < NX; i++)
				for (int j = 0; j < NY; j++) {
					Field[i][j].updateDxQ();
				}
		}
		else {
			for (int i = 0; i < NX; i++)
				for (int j = 0; j < NY; j++) {
					Field[i][j].updateDxEuler();
				}
		}
	}
	else {
		for (int i = 0; i < NX; i++)
			for (int j = 0; j < NY; j++) {
				Field[i][j].updateFastEuler(parity);
			}

	}
	// STOP TIMING
	auto tp2 = std::chrono::high_resolution_clock::now();
	elapsedTime = tp2 - tp1;

}

void SynchronizeFromFastToSlow() {
	for (int i = 0; i < NX; i++)
		for (int j = 0; j < NY; j++) {
			Field[i][j].dx[1] = Field[i][j].dx[0] - Field[i][j].v*dt;
			Field[i][j].clearQuad();

		}
}
void SynchronizeFromSlowToFast() {
	for (int i = 0; i < NX; i++)
		for (int j = 0; j < NY; j++) {
			Field[i][j].dx[0] = Field[i][j].dx[1];
		}
}

//void CreateSpeedGradient(int x, int y, int sizeX, int sizeY, double kmin, double kmax) {
//
//	for (int i = x; i < sizeX; i++)
//		for (int j = y; j < sizeY; j++) {
//			double X = map(j, 0, SIZEY, -2.1, 2.1);
//			O[i][j].k = map(-exp(-X * X), 0, -1, kmax, kmin);
//			if ((i == sizeX - 1 && j < sizeY - 50) || j == y || j == sizeY - 1)
//				(O[i][j]).active = true;
//		}
//}


double maxC = 0;

template <typename T>
void resetArray(std::array<std::array<T, NY>, NX> &X) {
	for (int i = 0; i < NX; i++)
		for (int j = 0; j < NY; j++) {
			//if (i < 2 || i > SIZEX - 2 || j < 2 || j > SIZEY - 2) {
			//	X[i][j] = T(
			//		0,
			//		(i+j)%2==0 ? 0 : 1);
			//	//if((i + j) % 2 == 0)
			//		X[i][j].active = false;
			//}
			//else {
				X[i][j] = T(
					(((j - NY / 2) * (j - NY / 2) + (i - NX*3 / 4) * (i - NX*3 / 4) < 60)) ? 4 : 0,
					isCircle(10, i - NX + 3 * NX / 4, j - NY / 2) ? 0.5 : 1);
			//}
		}
	//double kmax = 1;
	//double kmin = 0.4;
	//int x = 0;
	//int y = 0;
	//int sizeX = SIZEX * 2 / 3;
	//int sizeY = SIZEY ;
	//double alpha = 1;
	//for (int i = x; i < sizeX; i++)
	//	for (int j = y; j < sizeY; j++) {
	//		double Xc = map(j, 0, SIZEY, -2.1, 2.1);
	//		X[i][j].c = map(1-1./2*alpha*alpha*(Xc * Xc), 0,1, kmax, kmin);
	//		if ((i == sizeX - 1 && j < sizeY - 50) || j == y || j == sizeY - 1)
	//			(X[i][j]).active = true;
	//	}
	maxC = 1;
}




double maxDx(std::array<std::array<cell, NY>, NX>& X) {
	double max = 0;
	for (auto x : X)
		for (auto c : x)
			if (abs(c.dx[0]) > max)
				max = abs(c.dx[0]);
	return max;
}
bool AllZeros(std::array<std::array<cell, NY>, NX>& X) {
	for (auto x : X)
		for (auto c : x)
			if (abs(c.dx[0]) > 0.01)
				return false;
	return true;
}


class Example : public olc::PixelGameEngine
{
public:
	Example()
	{

		sAppName = "Cells";

	}

public:
	bool OnUserCreate() override
	{
		SetPixelMode(olc::Pixel::ALPHA);
		resetArray(Field);
		Join(Field);
		return true;
	}

	vd2d CheckPos = vd2d(0, 0);

	int parity = 0;
	bool OnUserUpdate(float fElapsedTime) override {
		Clear(olc::Pixel(0, 0, 0));
		int x = 0;
		int y = 0;

		for (const auto& vc : Field){
			for (const auto& C : vc) {
				olc::Pixel colour = olc::Pixel(enclose(abs(C.dx[parity]) * 255, 255), 0, enclose((C.dx[parity] + 1) * 255 / 2, 255), enclose(map(C.c, 0, 1, 0, 255), 255));
				FillRect(map(x, 0, NX, 0, SIZEX)-1, map(y, 0, NY, SIZEY, 0) - ((1.0 * SIZEY / NY) + 1), (1.0 * SIZEX / NX) + 1, ((1.0 * SIZEY / NY) + 1), colour);
				//Draw(x, y+SIZEY,  olc::Pixel(!CEuler[x][y].Q*255, CEuler[x][y].Q * 255,0));
				y++;
			}
		x++;
		y = 0;
		}
		Update(parity);
		if(!QuadraturePosition && !QuadratureVelocity)
			parity = 1 - parity;

		if (GetMouse(0).bPressed) {
			vd2d pos = GetMousePos();
		}
		if (GetMouse(1).bHeld) {
			vd2d pos = GetMousePos();
			double dx = Field[pos.x][pos.y].grad[0];
			double dy = Field[pos.x][pos.y].grad[1];
			//std::cout << dx << ", " << dy << "\n";
			DrawLine(pos, pos + vd2d(dx, dy) * 100);
		}
		if (GetKey(olc::Key::R).bPressed) {
			resetArray(Field);
			Join(Field);
		}

		if (GetKey(olc::Key::K1).bPressed) {
			QuadraturePosition = false;
			QuadratureVelocity = false;
			dt = h/1.5;
			SynchronizeFromFastToSlow();
		}
		if (GetKey(olc::Key::K2).bPressed) {
			QuadraturePosition = false;
			QuadratureVelocity = true;
			dt = h / 3.8;
			parity = 0;
			resetArray(Field);
			Join(Field);
			
		}
		if (GetKey(olc::Key::K3).bPressed) {
			QuadraturePosition = true;
			QuadratureVelocity = true;
			dt = h / 10;
			parity = 0;
			resetArray(Field);
			Join(Field);
		}

		DrawString(20, 20, std::to_string(elapsedTime.count()));
		DrawString(20,40, "dt "+std::to_string(dt));
		DrawString(20,60, "dx "+std::to_string(h));
		DrawString(20,80, "dx/dt"+std::to_string(h/dt));


		return true;
	}
};


int main(int argc, char* argv[])
{
	Example demo;
	if (demo.Construct(SIZEX, SIZEY, 1,1, false))
		demo.Start();
	return 0;


}