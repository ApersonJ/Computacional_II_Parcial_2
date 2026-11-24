#ifndef PARTICLE
#define PARTICLE

class Particle {
private:
    double x, v, m;

public:
    void Inicie(double x0, double v0, double mass){
		x = x0;
		v = v0;
		m = mass;
	}

    double momentum() const {
        return m * v;
    }

    double kineticEnergy() const {
        return 0.5 * m * v * v;
    }

    void setPosition(double xNew) {
        x = xNew;
    }

    void setVelocity(double vNew) {
        v = vNew;
    }

    double getX() const { return x; }
    double getV() const { return v; }
    double getM() const { return m; }
};

#endif
