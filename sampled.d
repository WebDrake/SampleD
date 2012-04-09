import std.random;
import std.stdio;

interface Sampler {
	size_t select();
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

			do {
				V = uniform(0.0,1.0);
			} while(V==0.0);

			while (quot > V) {
				++S;
				quot *= (top - S) / (recordsRemaining - S);
			}
		}

		return S;
	}
}

void main() {
	auto s = new VitterA(100,5);
	foreach(size_t i; 0..5) {
		writeln("Selected record: ", s.select());
	}
}
