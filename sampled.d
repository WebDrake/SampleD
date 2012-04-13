import std.math;
import std.random;
import std.stdio;
import std.c.time;


class SamplingAlgorithm(UniformRNG)
{
	this(size_t records, size_t sample, ref UniformRNG urng)
	{
		_recordsRemaining = _recordsTotal = records;
		_sampleRemaining = _sampleTotal = sample;
		_currentRecord = 0;
	}
	
	size_t select(ref UniformRNG urng)
	in
	{
		// No one should actually call this function.
		// Raise a big fuss if they do!!
		assert(false);
	}
	body
	{
		return 0;
	}

	final size_t recordsRemaining() { return _recordsRemaining; }
	final size_t recordsTotal() { return _recordsTotal; }

	final size_t sampleRemaining() { return _sampleRemaining; }
	final size_t sampleTotal() { return _sampleTotal; }

protected:
	size_t _currentRecord = 0;
	size_t _recordsRemaining;
	size_t _recordsTotal;
	size_t _sampleRemaining;
	size_t _sampleTotal;
}


class SamplingAlgorithmS(UniformRNG): SamplingAlgorithm!UniformRNG
{ 
	this(size_t records, size_t sample, ref UniformRNG urng)
	{
		super(records,sample,urng);
	}

	final size_t select(ref UniformRNG urng)
	in
	{
		assert(_sampleRemaining > 0);
		assert(_recordsRemaining >= _sampleRemaining);
	}
	body
	{
		//writeln("Starting S");
		while(true) {
			size_t r = uniform(0, _recordsRemaining, urng);

			if(r < _sampleRemaining) {
				size_t selectedRecord = _currentRecord;
				_sampleRemaining--;
				_recordsRemaining--;
				_currentRecord++;
				//writeln("\treturning ",selectedRecord,
				//        "\tremaining records: ", _recordsRemaining,
				//        "\tremaining to sample: ", _sampleRemaining,
				//        "\tnext record to consider: ", _currentRecord);
				return selectedRecord;
			} else {
				_recordsRemaining--;
				_currentRecord++;
			}
		}
	}
}


class SamplingAlgorithmSkip(UniformRNG): SamplingAlgorithm!UniformRNG
{
	this(size_t records, size_t sample, ref UniformRNG urng)
	{
		super(records,sample,urng);
	}
	
	final size_t select(ref UniformRNG urng)
	in
	{
		assert(_sampleRemaining > 0);
		assert(_recordsRemaining >= _sampleRemaining);
	}
	body
	{
		immutable size_t S = skip(urng);
		immutable size_t selectedRecord = _currentRecord + S;

		_sampleRemaining--;
		_recordsRemaining -= (S+1);
		_currentRecord += (S+1);
		
		return selectedRecord;
	}

protected:
	size_t skip(ref UniformRNG urng)
	{
		return 0;
	}
}


class SamplingAlgorithmA(UniformRNG): SamplingAlgorithmSkip!UniformRNG
{
	this(size_t records, size_t sample, ref UniformRNG urng)
	{
		super(records,sample,urng);
	}
	
protected:
	override size_t skip(ref UniformRNG urng) {
		size_t S;
		double V, quot, top;

		if(_sampleRemaining==1)
			S = uniform(0, _recordsRemaining, urng);
		else {
			S = 0;
			top = _recordsRemaining - _sampleRemaining;
			quot = top / _recordsRemaining;

			V = uniform!("()")(0.0,1.0, urng);

			while (quot > V) {
				++S;
				quot *= (top - S) / (_recordsRemaining - S);
			}
		}

		return S;
	}
}


class SamplingAlgorithmD(UniformRNG): SamplingAlgorithmA!UniformRNG
{
	this(size_t records, size_t sample, ref UniformRNG urng)
	{
		super(records,sample, urng);
		if( (alphaInverse * _sampleRemaining) > _recordsRemaining) {
			useVitterA = true;
		} else {
			Vprime = newVprime(_sampleRemaining, urng);
			useVitterA = false;
		}
	}
	
protected:
	override size_t skip(ref UniformRNG urng)
	{
		size_t S;
		size_t top, t, limit;
		size_t qu1 = 1 + _recordsRemaining - _sampleRemaining;
		double X, y1, y2, bottom;

		if(useVitterA)
			return super.skip(urng);
		else if ( (alphaInverse * _sampleRemaining) > _recordsRemaining) {
			useVitterA = true;
			return super.skip(urng);
		} else if ( _sampleRemaining > 1 ) {
			while(1) {
				for(X = _recordsRemaining * (1-Vprime), S = cast(size_t) trunc(X);
				    S >= qu1;
				    X = _recordsRemaining * (1-Vprime), S = cast(size_t) trunc(X)) {
					Vprime = newVprime(_sampleRemaining, urng);
				}

				y1 = (uniform!("()")(0.0,1.0, urng) * (cast(double) _recordsRemaining) / qu1) ^^ (1.0/(_sampleRemaining - 1));

				Vprime = y1 * ((-X/_recordsRemaining)+1.0) * ( qu1/( (cast(double) qu1) - S ) );

				if(Vprime > 1.0) {
					y2 = 1.0;
					top = _recordsRemaining -1;

					if(_sampleRemaining > (S+1) ) {
						bottom = _recordsRemaining - _sampleRemaining;
						limit = _recordsRemaining - S;
					} else {
						bottom = _recordsRemaining - (S+1);
						limit = qu1;
					}

					for( t=_recordsRemaining-1; t>=limit; --t) {
						y2 *= top/bottom;
						top--;
						bottom--;
					}

					if( (_recordsRemaining/(_recordsRemaining-X)) < (y1 * (y2 ^^ (1.0/(_sampleRemaining-1)))) ) {
						Vprime = newVprime(_sampleRemaining, urng);
					} else {
						Vprime = newVprime(_sampleRemaining-1, urng);
						return S;
					}
				} else {
					return S;
				}
			}
		} else {
			return cast(size_t) trunc(_recordsRemaining * Vprime);
		}
	}
private:
	immutable ushort alphaInverse = 13;
	double Vprime;
	bool useVitterA;
	double newVprime(size_t remaining, ref UniformRNG urng)
	{
		return ( uniform!("()")(0.0,1.0, urng) ^^ (1.0/remaining) );
	}
}


void sampling_test_simple(SamplerType, UniformRNG)
                         (size_t records, size_t samples, ref UniformRNG urng)
{
	auto s = new SamplerType(records,samples,urng);
	while(s.sampleRemaining > 0) {
		write("\trecord selected: ", s.select(urng), ".");
		write("\trecords remaining: ", s.recordsRemaining, ".");
		writeln("\tstill to sample: ", s.sampleRemaining, ".");
	}

	// Uncomment the following line if you want to force an assertion error ...
	// writeln("Let's sample 1 more for luck: ", s.select(urng), "\trecords remaining: ", s.recordsRemaining, "\tstill to sample: ", s.sampleRemaining);
}

void sampling_test_aggregate(SamplerType, UniformRNG)
                            (size_t records, size_t samples, ref UniformRNG urng,
                             size_t repeats=10000000, bool displayResults=true)
{
	double[] recordCount;
	clock_t start_time, end_time;

	writeln(SamplerType.classinfo.name, ", picking ", samples, " from ", records, ", ", repeats, " times.");

	recordCount.length = records;
	recordCount[] = 0.0;

	start_time = clock();

	foreach(size_t i; 0..repeats) {
		auto s = new SamplerType(records, samples, urng);

		while(s.sampleRemaining > 0) {
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
	sampling_test_simple!(SamplingAlgorithmA!Random,Random)(100,5,urng);
	writeln();
	
	writeln("Vitter's Algorithm D:");
	sampling_test_simple!(SamplingAlgorithmD!Random,Random)(100,5,urng);
	writeln();

	writeln("Fan et al./Jones' Algorithm S:");
	sampling_test_simple!(SamplingAlgorithmS!Random,Random)(100,5,urng);
	writeln();

	writeln("Now I'm going to again take samples of 5 from 100, but repeat the");
	writeln("process 10 million times.  This is basically a simple way of");
	writeln("checking for obvious bias -- we'll print out the number of times");
	writeln("each record is picked so you can see.");
	writeln();

	sampling_test_aggregate!(SamplingAlgorithmS!Random,Random)(100,5,urng,10000000,true);
	writeln();
	sampling_test_aggregate!(SamplingAlgorithmA!Random,Random)(100,5,urng,10000000,true);
	writeln();
	sampling_test_aggregate!(SamplingAlgorithmD!Random,Random)(100,5,urng,10000000,true);
	writeln();

	writeln("Now I'm going to take samples of 5 from 1000.  Notice how the time");
	writeln("required to carry out the sampling increases massively for Algorithm");
	writeln("A, but remains about the same for Algorithm D.");
	writeln();
	
	sampling_test_aggregate!(SamplingAlgorithmA!Random,Random)(1000,5,urng,10000000,false);
	writeln();
	
	sampling_test_aggregate!(SamplingAlgorithmD!Random,Random)(1000,5,urng,10000000,false);
	writeln();

	writeln("Now to finish up I'm going to sample 100,000 records from 10 million,");
	writeln("repeating the process 100 times.");
	writeln();
	sampling_test_aggregate!(SamplingAlgorithmA!Random,Random)(10000000,100000,urng,100,false);
	sampling_test_aggregate!(SamplingAlgorithmD!Random,Random)(10000000,100000,urng,100,false);
	writeln();

	writeln("That's all I'm demonstrating for now.  Hope you enjoyed the show!");
}
