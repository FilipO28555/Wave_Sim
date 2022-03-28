//#define OLC_PGE_APPLICATION
//#define _USE_MATH_DEFINES
//#include<Filip/olcPixelGameEngine.h>
//#include<Filip/LittleHelp.h>
//
//
//
//const int SIZEX = 600;
//const int SIZEY = 600;
//
//double dt = 0.1;  // timestep
//double pxpm = 1; //pixels per meter
//double h = 20;
//double speedC = 1; //speed of light in pixels per second
//
//using olc::vd2d;
//using std::cout;
//
//struct cell {
//	double dx,v,c; //wychylenie, prêdkoœæ wychylenia i prêdkoœæ œwiat³a
//	cell* N[4];   //s¹siedzi
//	double grad[2]; // = (gradx, grady). lewo - prawo, góra - dó³.
//	bool active = true;
//
//
//	//Quadratic extrapolation:
//	double a1, a2;
//	double v1, v2;
//	cell() {
//	}
//	cell(double DX, double C) {
//		a1 = 0;
//		a2 = 0;
//		v1 = 0;
//		v2 = 0;
//		dx = DX;
//		c = C;
//		v = 0;
//		for (int i = 0; i < 4; i++)
//			N[i] = 0;
//	}
//	cell(double DX, double C, double V) {
//		a1 = 0;
//		a2 = 0;
//		v1 = 0;
//		v2 = 0;
//		dx = DX;
//		c = C;
//		v = V;
//		for (int i = 0; i < 4; i++)
//			N[i] = 0;
//	}
//	//dt dependent on elapsed frame time. Using quadratic extrapolation assuming dt does not change.
//	void updateVQ(double dt) {
//		double DX = 0;
//		for (int i = 0; i < 4; i++)
//			if (N[i] && N[i]->active) {
//				double temp = (N[i]->dx - dx);
//				DX += temp;
//				grad[i / 2] = (1-2*(i%2))*temp;
//			}	
//			else {
//				double temp = v*h/c; //Dumpening on edges
//				DX -= temp;
//				grad[i / 2] = (1 - 2 * (i % 2)) * temp;
//			}
//		double a3 = (DX)/(h*h)*(c * c);
//		v += (3*a3-3*a2+a1) * dt;
//
//		a1 = a2;
//		a2 = a3;
//	}
//	//Quad also
//	void updateDxQ(double dt) {
//		dx += (3 * v - 3 * v2 + v1) * dt;
//		v1 = v2;
//		v2 = v;
//
//	}
//	void updateVEuler(double dt) {
//		double DX = 0;
//		for (int i = 0; i < 4; i++)
//			if (N[i] && N[i]->active)
//				DX += (N[i]->dx - dx);
//			else {
//				DX -= v * h; //Dumpening on edges
//			}
//		double a = c * c * (DX) / (h * h);
//		v += a * dt;
//
//	}
//	void updateDxEuler(double dt) {
//		dx += v*dt;
//
//	}
//
//};
//
//struct Tabcell {
//	double dx[2];
//	double c; //prêdkoœæ œwiat³a
//	Tabcell* N[4];   //s¹siedzi
//	bool active = true;
//
//
//	//Quadratic extrapolation:
//	double a1, a2;
//	double v1, v2;
//	Tabcell() {
//	}
//	Tabcell(double DX, double C) {
//		a1 = 0;
//		a2 = 0;
//		v1 = 0;
//		v2 = 0;
//		dx[0]=dx[1] = DX;
//		c = C;
//		
//		for (int i = 0; i < 4; i++)
//			N[i] = 0;
//	}
//	Tabcell(double DX, double C, double V) {
//		a1 = 0;
//		a2 = 0;
//		v1 = 0;
//		v2 = 0;
//		dx[0] = dx[1] = DX;
//		c = C;
//		for (int i = 0; i < 4; i++)
//			N[i] = 0;
//	}
//	//dt dependent on elapsed frame time. Using quadratic extrapolation assuming dt does not change.
//	void update(int par) {
//		double DX = 0;
//		for (int i = 0; i < 4; i++)
//			if (N[i] && N[i]->active)
//				DX += (N[i]->dx[par] - dx[par]);
//			else {
//				DX -= (dx[par]-dx[1-par])*dt*h * c; //Dumpening on edges
//			}
//		double a3 = (DX) / (h * h) * (c * c);
//		dx[1 - par] -= ((3 * a3 - 3 * a2 + a1)*dt) * dt + dx[1-par] ;
//
//		a1 = a2;
//		a2 = a3;
//	}
//};
//
//struct fastCell {
//	double dx[2];	    //
//	float c;		   //prêdkoœæ œwiat³a
//	fastCell* N[4];   //s¹siedzi
//	bool active = true;
//
//	fastCell() {
//	}
//	fastCell(double DX, double C) {
//		dx[0]=dx[1] = DX;
//		c = C;
//		for (int i = 0; i < 4; i++)
//			N[i] = 0;
//	}
//	fastCell(double DX, double C, double V) {
//		dx[0] = dx[1] = DX;
//		c = C;
//		for (int i = 0; i < 4; i++)
//			N[i] = 0;
//	}
//	void updateUsual(int par) {
//		double DX = 0;
//		for (int i = 0; i < 4; i++)
//			if (N[i] && N[i]->active)
//				DX += (N[i]->dx[par]);
//			else {
//				//DX -= (dx[par] - dx[1 - par]) / dt * h - dx[par]; //Dumpening on edges
//				DX -= (h/dt -1 )* dx[par] - dx[1 - par]*h/dt; //Dumpening on edges
//			}
//		double y2 = c * c * dt * dt / h / h;
//		dx[1 - par] = 2 * dx[par] - dx[1 - par] + (dt * dt * c * c / h / h) * (DX - 4 * dx[par]);
//
//	}
//};
//
//std::array<std::array<fastCell, SIZEY>, SIZEX> XEuler;
//std::array<std::array<cell, SIZEY>, SIZEX> XQuad;
//std::array<std::array<cell, SIZEY>, SIZEX> XDualQuad;
//
//template <typename T>
//void Join(std::array<std::array<T, SIZEY>, SIZEX> &X) {
//	for (int i = 0; i < SIZEX; i++)
//		for (int j = 0; j < SIZEY; j++) {
//			if (i + 1 < SIZEX)
//				X[i][j].N[0] = &X[i + 1][j]; //N[0] -> prawy
//			if (i - 1 >= 0)
//				X[i][j].N[1] = &X[i - 1][j]; //N[1] -> lewy
//			if (j + 1 < SIZEY)
//				X[i][j].N[2] = &X[i][j + 1]; //N[2] -> dolny
//			if (j - 1 >= 0)
//				X[i][j].N[3] = &X[i][j - 1]; //N[3] -> górny
//		}
//}
//
//std::chrono::duration<double> elapsedTimeEuler;
//std::chrono::duration<double> elapsedTimeQuad;
//std::chrono::duration<double> elapsedTimeDualQuad;
//
//void Update(double dt,int parity) {
//	// START TIMING
//	auto tp1 = std::chrono::high_resolution_clock::now();
//	//Quad - Euler
//	for (int i = 0; i < SIZEX; i++)
//		for (int j = 0; j < SIZEY; j++) {
//			XQuad[i][j].updateVQ(dt);
//		}
//	for (int i = 0; i < SIZEX; i++)
//		for (int j = 0; j < SIZEY; j++) {
//			XQuad[i][j].updateDxEuler(dt);
//
//		}
//	//for (int i = 0; i < SIZEX; i++)
//	//	for (int j = 0; j < SIZEY; j++) {
//	//		XQuad[i][j].update(parity);
//	//	}
//	// STOP TIMING
//	auto tp2 = std::chrono::high_resolution_clock::now();
//	elapsedTimeQuad = tp2 - tp1;
//
//
//
//
//	//DualQuad
//	tp1 = std::chrono::high_resolution_clock::now();
//	for (int i = 0; i < SIZEX; i++)
//		for (int j = 0; j < SIZEY; j++) {
//			XDualQuad[i][j].updateVQ(dt);
//		}
//	for (int i = 0; i < SIZEX; i++)
//		for (int j = 0; j < SIZEY; j++) {
//			XDualQuad[i][j].updateDxQ(dt);
//
//		}
//	tp2 = std::chrono::high_resolution_clock::now();
//	elapsedTimeDualQuad = tp2 - tp1;
//
//	////Euler - Euler
//	//tp1 = std::chrono::high_resolution_clock::now();
//	//if (parity) {
//	//	for (int i = 0; i < SIZEX; i++)
//	//		for (int j = 0; j < SIZEY; j++) {
//	//			XEuler[i][j].updateVEuler(dt);
//	//		}
//	//	for (int i = 0; i < SIZEX; i++)
//	//		for (int j = 0; j < SIZEY; j++) {
//	//			XEuler[i][j].updateDxEuler(dt);
//	//		}
//	//}
//	//tp2 = std::chrono::high_resolution_clock::now();
//	//elapsedTimeEuler = tp2 - tp1;
//
//	//Fast Euler
//	tp1 = std::chrono::high_resolution_clock::now();
//		for (int i = 0; i < SIZEX; i++)
//			for (int j = 0; j < SIZEY; j++) {
//				XEuler[i][j].updateUsual(parity);
//			}
//	tp2 = std::chrono::high_resolution_clock::now();
//	elapsedTimeEuler = tp2 - tp1;
//}
//
//
////void CreateSpeedGradient(int x, int y, int sizeX, int sizeY, double kmin, double kmax) {
////
////	for (int i = x; i < sizeX; i++)
////		for (int j = y; j < sizeY; j++) {
////			double X = map(j, 0, SIZEY, -2.1, 2.1);
////			O[i][j].k = map(-exp(-X * X), 0, -1, kmax, kmin);
////			if ((i == sizeX - 1 && j < sizeY - 50) || j == y || j == sizeY - 1)
////				(O[i][j]).active = true;
////		}
////}
//
//
//double maxC = 0;
//
//template <typename T>
//void resetArray(std::array<std::array<T, SIZEY>, SIZEX> &X) {
//	for (int i = 0; i < SIZEX; i++)
//		for (int j = 0; j < SIZEY; j++) {
//			//if (i < 2 || i > SIZEX - 2 || j < 2 || j > SIZEY - 2) {
//			//	X[i][j] = T(
//			//		0,
//			//		(i+j)%2==0 ? 0 : 1);
//			//	//if((i + j) % 2 == 0)
//			//		X[i][j].active = false;
//			//}
//			//else {
//				X[i][j] = T(
//					(((j - SIZEY / 2) * (j - SIZEY / 2) + (i - SIZEX*3 / 4) * (i - SIZEX*3 / 4) < 60)) ? 4 : 0,
//					isCircle(10, i - SIZEX + 3 * SIZEX / 4, j - SIZEY / 2) ? 0.5 : 1);
//			//}
//		}
//	//double kmax = 1;
//	//double kmin = 0.4;
//	//int x = 0;
//	//int y = 0;
//	//int sizeX = SIZEX * 2 / 3;
//	//int sizeY = SIZEY ;
//	//double alpha = 1;
//	//for (int i = x; i < sizeX; i++)
//	//	for (int j = y; j < sizeY; j++) {
//	//		double Xc = map(j, 0, SIZEY, -2.1, 2.1);
//	//		X[i][j].c = map(1-1./2*alpha*alpha*(Xc * Xc), 0,1, kmax, kmin);
//	//		if ((i == sizeX - 1 && j < sizeY - 50) || j == y || j == sizeY - 1)
//	//			(X[i][j]).active = true;
//	//	}
//	maxC = 1;
//}
//
//
//
////Convergence
//struct ConvInfo {
//public:
//	double h, dt;
//	bool Q = false;
//	ConvInfo() {
//		h = 0;
//		dt = 0;
//		Q = false;
//	}
//};
//std::array<std::array<ConvInfo, SIZEY>, SIZEX> CQuad;
//std::array<std::array<ConvInfo, SIZEY>, SIZEX> CEuler;
//std::array<std::array<ConvInfo, SIZEY>, SIZEX> CDualQuad;
//
//double maxDx(std::array<std::array<cell, SIZEY>, SIZEX>& X) {
//	double max = 0;
//	for (auto x : X)
//		for (auto c : x)
//			if (abs(c.dx) > max)
//				max = abs(c.dx);
//	return max;
//}
//bool AllZeros(std::array<std::array<cell, SIZEY>, SIZEX>& X) {
//	for (auto x : X)
//		for (auto c : x)
//			if (abs(c.dx) > 0.01)
//				return false;
//	return true;
//}
//void CheckConvergance(vd2d pos) {
//	resetArray(XQuad);
//	resetArray(XDualQuad);
//	resetArray(XEuler);
//	Join(XQuad);
//	Join(XDualQuad);
//	Join(XEuler);
//	dt = map(pos.x, 0, SIZEX, 0.01, 1);
//	h = map(pos.y, 0, SIZEY, 5,0.5);
//	int iter = 0;
//	while (++iter < 10) {
//		Update(dt,0);
//	}
//
////	CEuler[int(pos.x)][int(pos.y)].Q = (maxDx(XEuler) < 100)&& !AllZeros(XEuler);
////	CEuler[int(pos.x)][int(pos.y)].h = h;
////	CEuler[int(pos.x)][int(pos.y)].dt = dt;
////	CQuad[int(pos.x)][int(pos.y)].Q = (maxDx(XQuad) < 100)&& !AllZeros(XQuad);
////	CQuad[int(pos.x)][int(pos.y)].h = h;
////	CQuad[int(pos.x)][int(pos.y)].dt = dt;
////	CDualQuad[int(pos.x)][int(pos.y)].Q = (maxDx(XDualQuad) < 100)&& !AllZeros(XDualQuad);
////	CDualQuad[int(pos.x)][int(pos.y)].h = h;
////	CDualQuad[int(pos.x)][int(pos.y)].dt = dt;
//}
//void LogConvergence() {
//	std::ofstream LogE("ConverganceEuler.txt");
//	std::ofstream LogQ("ConverganceQuad.txt");
//	std::ofstream LogDQ("ConverganceDualQuad.txt");
//	for(int i=0;i<SIZEX;i++)
//		for (int j = 0; j < SIZEY; j++) {
//			//std::cout << i << " " << j << "\n";
//			LogE << CEuler[i][j].dt << ", " << CEuler[i][j].h << ", " << CEuler[i][j].Q << "\n";
//			LogQ << CQuad[i][j].dt << ", " << CQuad[i][j].h << ", " << CQuad[i][j].Q << "\n";
//			LogDQ << CDualQuad[i][j].dt << ", " << CDualQuad[i][j].h << ", " << CDualQuad[i][j].Q << "\n";
//	}
//}
//
//class Example : public olc::PixelGameEngine
//{
//public:
//	Example()
//	{
//
//		sAppName = "Cells";
//
//	}
//
//public:
//	bool OnUserCreate() override
//	{
//		SetPixelMode(olc::Pixel::ALPHA);
//		resetArray(XEuler);
//		resetArray(XQuad);
//		resetArray(XDualQuad);
//		Join(XEuler);
//		Join(XQuad);
//		Join(XDualQuad);
//		//LogConvergence();
//		return true;
//	}
//
//	bool finishedConvergence = true;
//	vd2d CheckPos = vd2d(0, 0);
//
//	int parity = 0;
//	bool OnUserUpdate(float fElapsedTime) override {
//		Clear(olc::Pixel(0, 0, 0));
//		int x = 0;
//		int y = 0;
//
//		for (const auto& vc : XEuler){
//			for (const auto& C : vc) {
//				Draw(x, y,  olc::Pixel(enclose(abs(C.dx[parity]) * 255,255),0, enclose((C.dx[parity] + 1) * 255 / 2,255), enclose(map(C.c, 0, 1, 0, 255),255)));
//				//Draw(x, y+SIZEY,  olc::Pixel(!CEuler[x][y].Q*255, CEuler[x][y].Q * 255,0));
//				y++;
//			}
//		x++;
//		y = 0;
//		}
//		 x = SIZEX;
//		 y = 0;
//		for (const auto& vc : XQuad) {
//			for (const auto& C : vc) {
//				Draw(x, y, olc::Pixel(enclose(abs(C.dx) * 255, 255), 0, enclose((C.dx + 1) * 255 / 2, 255), enclose(map(C.c, 0, 1, 0, 255), 255)));
//				//Draw(x, y + SIZEY, olc::Pixel(!CQuad[x-SIZEX][y].Q * 255, CQuad[x-SIZEX][y].Q * 255, 0));
//				y++;
//			}
//			x++;
//			y = 0;
//		}
//		x = 2*SIZEX;
//		y = 0;
//		for (const auto& vc : XDualQuad) {
//			for (const auto& C : vc) {
//				Draw(x, y, olc::Pixel(enclose(abs(C.dx) * 255, 255), 0, enclose((C.dx + 1) * 255 / 2, 255), enclose(map(C.c, 0, 1, 0, 255), 255)));
//				if (!((x+10) % 20) && !((y+10) % 20)) {
//					double dx = C.grad[0];
//					double dy = C.grad[1];
//					DrawLine(x, y,enclose(x+dx*50,3*SIZEX),enclose(y+dy*50,3*SIZEY), olc::Pixel(255, 255, 255, 200));
//				}
//				//Draw(x, y + SIZEY, olc::Pixel(!CDualQuad[x - 2*SIZEX][y].Q * 255, CDualQuad[x - 2*SIZEX][y].Q * 255, 0));
//				y++;
//			}
//			x++;
//			y = 0;
//		}
//		Update(dt, parity);// *fElapsedTime);
//		parity = 1 - parity;
//
//		if (!finishedConvergence) {
//			CheckConvergance(CheckPos);
//
//			CheckPos.x += 1;
//			if (CheckPos.x >= SIZEX) {
//				CheckPos.x = 0;
//				CheckPos.y++;
//				if (CheckPos.y > SIZEY) {
//					finishedConvergence = true;
//					LogConvergence();
//				}
//			}
//
//		}
//		if (GetMouse(0).bPressed) {
//			vd2d pos = GetMousePos();
//			CheckConvergance(pos);
//		}
//		if (GetMouse(1).bHeld) {
//			vd2d pos = GetMousePos();
//			double dx = XQuad[pos.x][pos.y].grad[0];
//			double dy = XQuad[pos.x][pos.y].grad[1];
//			//std::cout << dx << ", " << dy << "\n";
//			DrawLine(pos, pos + vd2d(dx, dy) * 100);
//		}
//		if (GetKey(olc::Key::R).bPressed) {
//			resetArray(XQuad);
//			resetArray(XDualQuad);
//			resetArray(XEuler);
//			Join(XQuad);
//			Join(XDualQuad);
//			Join(XEuler);
//		}
//
//
//		DrawString(20, 20, std::to_string(elapsedTimeEuler.count()));
//		DrawString(SIZEX+20, 20, std::to_string(elapsedTimeQuad.count()));
//		DrawString(2*SIZEX+20, 20, std::to_string(elapsedTimeDualQuad.count()));
//		DrawString(20,40, "dt "+std::to_string(dt));
//		DrawString(20,60, "dx "+std::to_string(h));
//		DrawString(20,80, "dx/dt"+std::to_string(h/dt));
//
//
//		return true;
//	}
//};
//
//
//int main(int argc, char* argv[])
//{
//	Example demo;
//	if (demo.Construct(3*SIZEX, SIZEY, 1,1, false))
//		demo.Start();
//	return 0;
//
//
//}