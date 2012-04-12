import std.math;
import std.random;
import std.stdio;
import std.c.time;

interface Sampler(UniformRandomNumberGenerator) {
	size_t select(ref UniformRandomNumberGenerator urng);
}


class Skipper(UniformRandomNumberGenerator) : Sampler!(UniformRandomNumberGenerator) {
	this(size_t records,
	     size_t sample,
	     ref UniformRandomNumberGenerator urng)
	{
		recordsRemaining = recordsTotal = records;
		sampleRemaining = sampleTotal = sample;
	}
	
	final size_t select(ref UniformRandomNumberGenerator urng)
	{
		if(recordsRemaining < 1) {
			// Throw error: no more records to sample.
		} else if (sampleRemaining < 1) {
			// Throw error: entire sample already aquired.
		}

		size_t S = skip(urng);

		sampleRemaining--;
		recordsRemaining -= (S+1);
		currentRecord += S;
		
		return currentRecord++;
	}

	final size_t records_remaining() { return recordsRemaining; }
	final size_t records_total() { return recordsTotal; }

	final size_t sample_remaining() { return sampleRemaining; }
	final size_t sample_total() { return sampleTotal; }

protected:
	size_t currentRecord = 0;
	size_t recordsRemaining;
	size_t recordsTotal;
	size_t sampleRemaining;
	size_t sampleTotal;
	size_t skip(ref UniformRandomNumberGenerator urng)
	{
		return 0;
	}
}

class VitterA(UniformRandomNumberGenerator) : Skipper!UniformRandomNumberGenerator {
	this(size_t records,
	     size_t sample,
	     ref UniformRandomNumberGenerator urng)
	{
		super(records,sample,urng);
	}
protected:
	override size_t skip(ref UniformRandomNumberGenerator urng) {
		size_t S;
		double V, quot, top;

		if(sampleRemaining==1)
			S = uniform(0, recordsRemaining, urng);
		else {
			S = 0;
			top = recordsRemaining - sampleRemaining;
			quot = top / recordsRemaining;

			V = uniform!("()")(0.0,1.0, urng);

			while (quot > V) {
				++S;
				quot *= (top - S) / (recordsRemaining - S);
			}
		}

		return S;
	}
}

class VitterD(UniformRandomNumberGenerator) : VitterA!(UniformRandomNumberGenerator) {
	this(size_t records,
	     size_t sample,
	     ref UniformRandomNumberGenerator urng)
	{
		super(records,sample, urng);
		if( (alphaInverse * sampleRemaining) > recordsRemaining) {
			useVitterA = true;
		} else {
			Vprime = newVprime(sampleRemaining, urng);
			useVitterA = false;
		}
	}
	
protected:
	override size_t skip(ref UniformRandomNumberGenerator urng) {
		size_t S;
		size_t top, t, limit;
		size_t qu1 = 1 + recordsRemaining - sampleRemaining;
		double X, y1, y2, bottom;

		if(useVitterA)
			return super.skip(urng);
		else if ( (alphaInverse * sampleRemaining) > recordsRemaining) {
			useVitterA = true;
			return super.skip(urng);
		} else if ( sampleRemaining > 1 ) {
			while(1) {
				for(X = recordsRemaining * (1-Vprime), S = cast(size_t) trunc(X);
				    S >= qu1;
				    X = recordsRemaining * (1-Vprime), S = cast(size_t) trunc(X)) {
					Vprime = newVprime(sampleRemaining, urng);
				}

				y1 = pow( (uniform!("()")(0.0,1.0, urng) * (cast(double) recordsRemaining) / qu1),
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
						Vprime = newVprime(sampleRemaining, urng);
					} else {
						Vprime = newVprime(sampleRemaining-1, urng);
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
	double newVprime(size_t remaining,
	                 ref UniformRandomNumberGenerator urng)
	{
		return pow ( uniform!("()")(0.0,1.0, urng), 1.0/remaining );
	}
}


void sampling_test_simple(SamplerType, UniformRandomNumberGenerator)
(size_t records, size_t samples, ref UniformRandomNumberGenerator urng)
{
	auto s = new SamplerType(records,samples,urng);
	while(s.sample_remaining() > 0) {
		write("\trecord selected: ", s.select(urng), ".");
		write("\trecords remaining: ", s.records_remaining(), ".");
		writeln("\tstill to sample: ", s.sample_remaining(), ".");
	}
}

void sampling_test_aggregate(SamplerType, UniformRandomNumberGenerator)
(size_t records, size_t samples, ref UniformRandomNumberGenerator urng, size_t repeats=10000000, bool displayResults=true)
{
	double[] recordCount;
	clock_t start_time, end_time;

	writeln(SamplerType.classinfo.name, ", picking ", samples, " from ", records, ", ", repeats, " times.");

	recordCount.length = records;
	recordCount[] = 0.0;

	start_time = clock();

	foreach(size_t i; 0..repeats) {
		auto s = new SamplerType(records, samples, urng);

		while(s.sample_remaining() > 0) {
			recordCount[s.select(urng)]++;
		}
	}

	end_time = clock();

	if(displayResults) {
		foreach(size_t record, double count; recordCount)
			writeln("\trecord ", record, " was picked ", count, " times.");
	}

	writeln("\t\tSampling completed in ",(cast(double) (end_time - start_time))/CLOCKS_PER_SEC, " seconds.");
}

void main() {
	auto urng = Random(23);

	writeln("Hello!  I'm a very simple and stupid implementation of some clever");
	writeln("sampling algorithms.");
	writeln();
	writeln("What I'm going to do first is just take a sample of 5 records out of");
	writeln("100 using each of the algorithms.");
	writeln();
	
	writeln("Vitter's Algorithm A:");
	sampling_test_simple!(VitterA!Random,Random)(100,5,urng);
	writeln();
	
	writeln("Vitter's Algorithm D:");
	sampling_test_simple!(VitterD!Random,Random)(100,5,urng);
	writeln();

	writeln("Now I'm going to again take samples of 5 from 100, but repeat the");
	writeln("process 10 million times.  This is basically a simple way of");
	writeln("checking for obvious bias -- we'll print out the number of times");
	writeln("each record is picked so you can see.");
	writeln();

	sampling_test_aggregate!(VitterA!Random,Random)(100,5,urng,10000000,true);
	writeln();
	sampling_test_aggregate!(VitterD!Random,Random)(100,5,urng,10000000,true);
	writeln();

	writeln("Now I'm going to take samples of 5 from 1000.  Notice how the time");
	writeln("required to carry out the sampling increases massively for Algorithm");
	writeln("A, but remains about the same for Algorithm D.");
	writeln();
	
	sampling_test_aggregate!(VitterA!Random,Random)(1000,5,urng,10000000,false);
	writeln();
	
	sampling_test_aggregate!(VitterD!Random,Random)(1000,5,urng,10000000,false);
	writeln();

	writeln("Now to finish up I'm going to sample 100,000 records from 10 million,");
	writeln("repeating the process 100 times.");
	writeln();
	sampling_test_aggregate!(VitterA!Random,Random)(10000000,100000,urng,100,false);
	sampling_test_aggregate!(VitterD!Random,Random)(10000000,100000,urng,100,false);
	writeln();

	writeln("That's all I'm demonstrating for now.  Hope you enjoyed the show!");
}
