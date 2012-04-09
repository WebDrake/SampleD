import std.math;
import std.random;
import std.stdio;

interface Sampler {
	size_t select();
}

double uniform_pos(double a, double b) {
	double x;

	do {
		x = uniform(a,b);
	} while (x==0.0);

	return x;
}

class Skipper : Sampler {
	this(size_t records,
	     size_t sample) {
		recordsRemaining = recordsTotal = records;
		sampleRemaining = sampleTotal = sample;
	}
	
	final size_t select() {
		if(recordsRemaining < 1) {
			// Throw error: no more records to sample.
		} else if (sampleRemaining < 1) {
			// Throw error: entire sample already aquired.
		}

		size_t S = skip() + 1;

		sampleRemaining--;
		recordsRemaining -= S;

		currentRecord += S;
		return currentRecord;
	}
protected:
	size_t currentRecord = 0;
	size_t recordsRemaining;
	size_t recordsTotal;
	size_t sampleRemaining;
	size_t sampleTotal;
	size_t skip() {
		return 0;
	}
}

class VitterA : Skipper {
	this(size_t records,
	     size_t sample) {
		super(records,sample);
	}
protected:
	override size_t skip() {
		size_t S;
		double V, quot, top;

		if(sampleRemaining==1)
			S = uniform(0, recordsRemaining);
		else {
			S = 0;
			top = recordsRemaining - sampleRemaining;
			quot = top / recordsRemaining;

			V = uniform_pos(0.0,1.0);

			while (quot > V) {
				++S;
				quot *= (top - S) / (recordsRemaining - S);
			}
		}

		return S;
	}
}

class VitterD : VitterA {
	this(size_t records,
	     size_t sample) {
		super(records,sample);
		if( (alphaInverse * sampleRemaining) > recordsRemaining) {
			useVitterA = true;
		} else {
			Vprime = newVprime(sampleRemaining);
			useVitterA = false;
		}
	}
	
protected:
	override size_t skip() {
		size_t S;
		size_t top, t, limit;
		size_t qu1 = 1 + recordsRemaining - sampleRemaining;
		double X, y1, y2, bottom;

		if(useVitterA)
			return super.skip();
		else if ( (alphaInverse * sampleRemaining) > recordsRemaining) {
			useVitterA = true;
			return super.skip();
		} else if ( sampleRemaining > 1 ) {
			while(1) {
				for(X = recordsRemaining * (1-Vprime), S = cast(size_t) trunc(X);
				    S >= qu1;
				    X = recordsRemaining * (1-Vprime), S = cast(size_t) trunc(X)) {
					Vprime = newVprime(sampleRemaining);
				}

				y1 = pow( (uniform_pos(0.0,1.0) * (cast(double) recordsRemaining) / qu1),
				          (1.0/(sampleRemaining - 1)) );

				Vprime = y1 * ((-X/recordsRemaining)+1.0) * ( qu1/( (cast(double) qu1) - S ) );

				if(Vprime > 1.0) {
					y2 = 1.0;
					top = recordsRemaining -1;

					if(sampleRemaining > (S+1) ) {
						bottom = recordsRemaining - sampleRemaining;
						limit = recordsRemaining - S;
					} else {
						bottom = recordsRemaining - (S+1);
						limit = qu1;
					}

					for( t=recordsRemaining-1; t>=limit; --t)
						y2 *= top--/bottom--;

					if( (recordsRemaining/(recordsRemaining-X)) < (y1*pow(y2, 1.0/(sampleRemaining-1))) ) {
						Vprime = newVprime(sampleRemaining);
					} else {
						Vprime = newVprime(sampleRemaining-1);
						return S;
					}
				} else {
					return S;
				}
			}
		} else {
			return cast(size_t) trunc(recordsRemaining * Vprime);
		}
	}
private:
	immutable ushort alphaInverse = 13;
	double Vprime;
	bool useVitterA;
	double newVprime(size_t remaining) {
		return pow ( uniform_pos(0.0,1.0), 1.0/remaining );
	}
}

void main() {
	auto s = new VitterD(100,5);
	foreach(size_t i; 0..5) {
		writeln("Selected record: ", s.select());
	}
}
