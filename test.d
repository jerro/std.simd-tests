import std.stdio, std.traits, std.math, std.conv, std.typetuple, std.string,
    std.bigint, std.random, std.algorithm, std.range, std.getopt;
import std.simd;
import core.simd;

alias TypeTuple tt;

template isVector(T)
{
    static if(is(T == double2) || is(T == float4) ||
            is(T == long2) || is(T == ulong2) ||
            is(T == int4) || is(T == uint4) ||
            is(T == short8) || is(T == ushort8) ||
            is(T == byte16) || is(T == ubyte16) || is(T == void16))
        enum bool isVector = true;
    else
        enum bool isVector = false;
}

template PromotionOf(T)
{
    static if(is(T == int4))
        alias long2 PromotionOf;
    else static if(is(T == uint4))
        alias ulong2 PromotionOf;
    else static if(is(T == short8))
        alias int4 PromotionOf;
    else static if(is(T == ushort8))
        alias uint4 PromotionOf;
    else static if(is(T == byte16))
        alias short8 PromotionOf;
    else static if(is(T == ubyte16))
        alias ushort8 PromotionOf;
    else static if(is(T == int))
        alias long PromotionOf;
    else static if(is(T == uint))
        alias ulong PromotionOf;
    else static if(is(T == short))
        alias int PromotionOf;
    else static if(is(T == ushort))
        alias uint PromotionOf;
    else static if(is(T == byte))
        alias short PromotionOf;
    else static if(is(T == ubyte))
        alias ushort PromotionOf;
    else
        static assert(0, "Incorrect type");
}

template BaseType(V)
{
    alias typeof(V.array[0]) BaseType;
}

template nElements(V)
{
    enum nElements = V.sizeof / BaseType!(V).sizeof;
}

template staticIota(int start, int end, int stride = 1)
{
    static if(start >= end)
        alias tt!() staticIota;
    else
        alias tt!(start, staticIota!(start + stride, end, stride)) 
            staticIota;
}

template staticRepeat(int n, a...) if(a.length == 1)
{
    static if(n <= 0)
        alias tt!() staticRepeat;
    else
        alias tt!(a, staticRepeat!(n - 1, a)) staticRepeat;
}

auto vector(T, A...)(A a)
{
    alias Vector!(T[A.length]) V;
    V r;
    foreach(i,_; A)
        (cast(T*) &r)[i] = cast(T) a[i];
    
    return r;
}

auto repeated(V)(BaseType!V a)
{
    V r;
    r = a;
    return r;
}

auto eq(bool approx = false, T)(T a, T b) if(isIntegral!T || is(T == bool))
{
    return a == b;
}

auto eq(bool approx = false, T)(T a, T b) if(isFloatingPoint!T)
{
    if(a.isnan && b.isnan)
        return true;

    return feqrel(a, b) + 3 >= (approx ? T.mant_dig / 2 : T.mant_dig);
}

auto eq(bool approx = false, V)(V a, V b) if(isVector!V)
{
    alias typeof(b.array[0]) T;

    foreach(i; staticIota!(0, V.init.length))
        if(!eq!approx(a.array[i], b.array[i]))
            return false;

    return true;
}

/*auto eq(bool approx = false, V, V1)(V a, V1 bb) 
    if(isVector!V && isImplicitlyConvertible!(V1, V) && !is(V == V1))
{
    V b = bb;
    return eq(a, b);
}*/

auto eq(bool approx = false, V, V1)(V1 aa, V b) 
    if(isVector!V && isImplicitlyConvertible!(V1, V) && !is(V == V1))
{
    V a = aa;
    return eq(a, a);
}

template getString(a...)
{
    static if(a.length == 0)
        enum getString = "";
    else
    {
        static if(is(typeof(a[0]) == string))
            enum getString = `"` ~ a[0] ~ `" ` ~ getString!(a[1 .. $]);
        else
            enum getString = a[0].stringof ~ " " ~ getString!(a[1 .. $]);
    }
}

void print(T)(File f, T a)
{
    static if(is(typeof(a.array)))
        f.writeln(a.array);
    else
        f.writeln(a);
}

template group(a...){ alias a members; }

alias tt!(byte, ubyte, short, ushort, int, uint, long, ulong) integral;
alias tt!(float, double) floatingPoint;

bool listSymbols = false;

template test(alias op, bool approx = false, templateParams...)
{
    static if (templateParams.length == 0)
        alias tt!(SIMDVer.SSE42) tParams;
    else
        alias templateParams tParams;

    void test(A...)(A a)
    {
        alias a[0 .. $ - 1] params;
        alias a[$ - 1] correct;
        
        static if(is(typeof(op!tParams(params))))
        {
            if(listSymbols)       
            {
                //try to find the mangled name and print it
                static if(is(typeof(op!tParams)))
                {
                    alias op!(tParams) fun;
                    writeln(fun.mangleof);
                    return;
                }
                foreach(p; params)
                    static if(
                        is(typeof(op!(tParams, typeof(p)))) && 
                        is(typeof(op!(tParams, typeof(p))(params)) == 
                        typeof(op!tParams(params))))
                    {
                        alias op!(tParams, typeof(p)) fun;
                        writeln(fun.mangleof);
                        return;
                    }
            }
            else
            {
                auto r = op!tParams(params);   

                if(!eq!approx(r, correct))
                {
                    stderr.writeln("Template parameters:");
                    foreach(p; tParams)
                        static if(is(typeof(p)))
                            stderr.print(p);

                    stderr.writeln("Parameters:");
                    foreach(p; params)
                        stderr.print(p);
                    
                    stderr.writeln("Parameter types:");
                    foreach(p; params)
                        stderr.print(typeof(p).stringof);
                    
                    stderr.writeln("Result: ");
                    stderr.print(r);
                    
                    stderr.writeln("Correct result: ");
                    stderr.print(correct);
     
                    assert(false, format(
                        "Function %s returned an incorrect result", 
                        op.stringof));
                }
            }           
        }
        else
            pragma(msg, 
                "Failed to compile: " ~ op.stringof ~
                "   params: " ~ typeof(params).stringof ~ 
                "   tParams:  " ~ getString!tParams);
    }
}

void testElementWise(
    alias finitGroup, alias opsGroup, bool approx = false, 
    SIMDVer ver = SIMDVer.SSE42, int vecBytes = 16, Scal...)()
{
    static if(Scal.length == 0)
        alias tt!(floatingPoint, integral) Scalars;
    else 
        alias Scal Scalars; 
   
    alias finitGroup.members finit;
    alias opsGroup.members ops;
 
    enum nParams = finit.length;
    enum nOps = ops.length / 2;
 
    foreach(T; Scalars)
    {
        enum n = vecBytes / T.sizeof;
        alias Vector!(T[n]) V;
        staticRepeat!(nParams, V) params;

        foreach(i, _; params)
            foreach(long j; 0 .. n)
                (cast(T*) &params[i])[j] = cast(T) finit[i](j);

        foreach(i; staticIota!(0, 2 * nOps, 2) )
        {
            V correct;
            foreach(j; staticIota!(0, n))
            {
                staticRepeat!(nParams, T) scalarParams;
                foreach(k, _2; scalarParams)
                    scalarParams[k] = params[k].array[j];
               
                (cast(T*) &correct)[j] = cast(T)ops[i + 1](scalarParams);
            }
           
            alias ops[i] op; 

            test!(op, approx, ver)(params, correct);
        }
    }
}

template opSaturate(string op)
{
    auto opSaturate(T)(T a, T b)
    {
        static if(isIntegral!T)
        {
            auto ba = BigInt(a);
            auto bb = BigInt(b);
            auto br = mixin("ba" ~ op ~ "bb");
            return 
                br > T.max ? T.max :
                br < T.min ? T.min : cast(T) br.toLong();
        }
        else
            return a + b;
    }
}

@property intBitcast(T)(T a)
{
    static if(is(T == float))
        return *cast(uint*)& a;
    else static if(is(T == double))
        return *cast(ulong*)& a;
    else static if(isIntegral!T)
        return a;
    else
        static assert(0); 
}

@property bitcast(R, T)(T a)
{
    return *cast(R*)& a;
}

auto randomSwizzleString(int n, int seed)
{
    assert(n <= 16);
    auto rng = Xorshift128(seed);
    char[16] r;
    foreach(ref e; r)
        e = "0123456789ABCDEF"[uniform(0, n, rng)];

    return r[0 .. n].idup;
}

template swizzleStrings(int nElements, int maxNStrings)
{
    auto f() 
    {
        assert(nElements <= 16);
        bool[string] table;
        foreach(i; 0 .. maxNStrings)
            table[`"` ~ randomSwizzleString(nElements, i) ~ `"`] = true;

        return "alias tt!(" ~ table.keys.join(", ").idup ~ ") swizzleStrings;";
    }
    
    mixin(f);
}

@property ptr_(V)(ref V v)
{
    return cast(BaseType!(V)*) &v;
}

auto simpleSwizzle(V)(string ind, V v)
{
    V r;
    foreach(i, char c; ind)
        r.ptr_[i] = v.array["0123456789ABCDEF".countUntil(c)]; 
    
    return r;
}

auto testStoreScalar(SIMDVer ver, V)(V vec)
{
    BaseType!V s;
    storeScalar!ver(vec, &s);
    return s;
}

auto testStoreUnaligned(SIMDVer ver, V)(V v)
{
    BaseType!(V)[nElements!V + 1] mem;
    storeUnaligned!ver(v, &mem[1]);
    V r;
    r.ptr_[0 .. nElements!V] = mem[1 .. $];
    return r;
}

void main(string[] args)
{
    getopt(args, "l", &listSymbols); 

    // precise unary floating point operations:
    testElementWise!(
        group!(a => 0.55 * a + 0.51), 
        group!(
            std.simd.floor, std.algorithm.floor,
            std.simd.ceil, std.algorithm.ceil,
            std.simd.round, std.algorithm.round,
            std.simd.trunc, std.algorithm.trunc),
        false, SIMDVer.SSE42, 16, floatingPoint)();

    // imprecise unary floating point operations:
    testElementWise!(
        group!(a => 0.55 * a + 0.51), 
        group!(
            rcp, a => cast(typeof(a)) 1 / a,
            rcpEst, a => cast(typeof(a)) 1 / a,
            rsqrt, a => cast(typeof(a)) 1 / std.math.sqrt(a),
            std.simd.sqrt, std.math.sqrt,
            rsqrtEst, a => cast(typeof(a)) 1 / std.math.sqrt(a),
            std.simd.sqrtEst, std.math.sqrt),
        true, SIMDVer.SSE42, 16, floatingPoint)();

    // unary operations
    testElementWise!(
        group!(a => 2 * a - 1), 
        group!(
            comp, a => (~a.intBitcast).bitcast!(typeof(a)),
            neg, a => -a,
            std.simd.abs, a => a < 0 ? -a : a))();

    // binary operations:
    testElementWise!(
        group!(a => a + 2, a => 2 * a + 1), 
        group!(
            add, (a, b) => a + b,
            sub, (a, b) => a - b,
            mul, (a, b) => a * b,
            div, (a, b) => a / b,
            divEst, (a, b) => a / b,
            std.simd.min, std.algorithm.min,
            std.simd.max, std.algorithm.max,
            or, (a, b) => (a.intBitcast | b.intBitcast).bitcast!(typeof(a)),
            nor, (a, b) => (~(a.intBitcast | b.intBitcast)).bitcast!(typeof(a)),
            and, (a, b) => (a.intBitcast & b.intBitcast).bitcast!(typeof(a)),
            nand, (a, b) => (~(a.intBitcast & b.intBitcast)).bitcast!(typeof(a)),
            andNot, (a, b) => (a.intBitcast & ~b.intBitcast).bitcast!(typeof(a)),
            xor, (a, b) => (a.intBitcast ^ b.intBitcast).bitcast!(typeof(a))))();
 
    // integral binary operations (onlys supported for byte16, ubyte16,
    // short8 and ushort8 on x86):
    testElementWise!(
        group!(a => a + 2, a => 2 * a + 1), 
        group!(
            addSaturate, opSaturate!"+",
            subSaturate, opSaturate!"-"),
        false, SIMDVer.SSE42, 16, tt!(byte, ubyte, short, ushort))();

    // three operand operations:
    testElementWise!(
        group!(a => a + 2, a => 2 * a + 1, a => 3 - a),
        group!(
            madd, (a, b, c) => a * b + c,
            msub, (a, b, c) => a * b - c,
            nmadd, (a, b, c) => - a * b + c,
            nmsub, (a, b, c) => - a * b - c,
            lerp, (a, b, c) => a + c * (b - a),
            clamp, (a, b, c) => std.algorithm.max(a, std.algorithm.min(b, c))))();

    // four operand operations:
    testElementWise!(
        group!(a => a + 2, a => 2 * a + 1, a => 3 - a, a => 2),
        group!(
            selectEqual, (a, b, c, d) => a == b ? c : d,
            selectNotEqual, (a, b, c, d) => a != b ? c : d,
            selectGreater, (a, b, c, d) => a > b ? c : d,
            selectGreaterEqual, (a, b, c, d) => a >= b ? c : d,
            selectLess, (a, b, c, d) => a < b ? c : d,
            selectLessEqual, (a, b, c, d) => a <= b ? c : d))();

    enum vecBytes = 16;

    foreach(T; tt!(floatingPoint, integral))
    {
        enum n = vecBytes / T.sizeof; 
        alias Vector!(T[n]) V;       
        alias staticIota!(1, n + 1) seq;
        
        Vector!(ubyte[V.sizeof]) mask;
        V v, v1, v2, v3; 
        V ones = 1;
        V zeros = 0;
        V correct;
        V vseq = vector!T(seq);

        T s = 1;
        correct = 0;
        correct.ptr_[0] = 1;
        test!(loadScalar, false, V, SIMDVer.SSE42)(&s, correct);
       
        T[n + 1] a = [0, seq];
        test!(loadUnaligned, false, V, SIMDVer.SSE42)(&a[1], vector!T(seq));

        test!getScalar(vector!T(seq), to!T(1));
        test!testStoreScalar(vector!T(seq), to!T(1));
        test!testStoreUnaligned(vector!T(seq), vector!T(seq));

        correct = 0;
        correct.ptr_[0] = 1;
        test!setX(zeros, ones, correct);
        test!getX(correct, ones);

        static if(n >= 2)
        {
            correct = 0;
            correct.ptr_[1] = 1;
            test!setY(zeros, ones, correct);
            test!getY(correct, ones);
        }
        
        static if(n >= 3)
        {
            correct = 0;
            correct.ptr_[2] = 1;
            test!setZ(zeros, ones, correct);
            test!getZ(correct, ones);
        }

        static if(n >= 4)
        {
            correct = 0;
            correct.ptr_[3] = 1;
            test!setW(zeros, ones, correct);
            test!getW(correct, ones);
        }
        
        foreach(ind; swizzleStrings!(n, 50))
            test!(swizzle, false, ind, SIMDVer.SSE42)(
                vseq, simpleSwizzle(ind, vseq));
        
        // broadcasts
        foreach(i; staticIota!(0, n))
        {
            enum swiz = "" ~ "0123456789abcdef"[i];
            v = 0;
            v.ptr_[i] = 1;
            test!(swizzle, false, swiz, SIMDVer.SSE42)(v, ones); 
        }

        {
            alias staticIota!(0, 2 * n) both;
            v1 = vector!T(both[0 .. n]);
            v2 = vector!T(both[n .. 2 * n]);
            auto result = new T[2 * n];
            foreach(i; both)
               result[i] =  n * (i & 1) + (i >> 1);
            correct.ptr_[0 .. n] = result[0 .. n];
            test!interleaveLow(v1, v2, correct);
            correct.ptr_[0 .. n] = result[n .. $];
            test!interleaveHigh(v1, v2, correct);
        }

        static if(is(PromotionOf!T))
        {
            alias PromotionOf!T P;
            alias Vector!(P[n / 2]) VP;
            VP vp1, vp2;

            v = vector!T(staticIota!(0, n));
            vp1 = vector!P(staticIota!(0, n / 2));
            vp2 = vector!P(staticIota!(n / 2, n));
            test!unpackLow(v, vp1); 
            test!unpackHigh(v, vp2);
            test!pack(vp1, vp2, v);
            // Does not make sense for floating point and only works with
            // vectors with 16 bit or 8 bit elements on x86.
            static if(isIntegral!T && T.sizeof <= 2)
                test!packSaturate(vp1, vp2, v); // TODO: check that it saturates
        }

        static if(n == 4)
        {
            v = vector!T(seq);
            test!toFloat(v, vector!float(seq));
            test!toInt(v, vector!int(seq));
        }

        static if(n == 2)
            test!toDouble(vector!T(seq), vector!double(seq));
    
        static if(isIntegral!T && T.sizeof > 1)
        {
            v1 = vector!T(staticIota!(0, 2 * n, 2));
            v2 = 0;
            v2.ptr_[0] = 1;

            test!shiftLeft(v1, v2, vector!T(staticIota!(0, 4 * n, 4))); 
            test!(shiftLeftImmediate, false, 1, SIMDVer.SSE42)(
                v1, vector!T(staticIota!(0, 4 * n, 4))); 
            
            test!shiftRight(v1, v2, vector!T(staticIota!(0, n))); 
            test!(shiftRightImmediate, false, 1, SIMDVer.SSE42)(
                v1, vector!T(staticIota!(0, n)));
        } 

        mask = 0;
        foreach(i; 0 .. T.sizeof)
            mask.ptr_[i] = ubyte.max;
        v1 = vector!T(staticIota!(50, 50 + n));
        v2 = vector!T(seq);
        correct = v2;
        correct.ptr_[0] = v1.array[0];
        test!(std.simd.select)(mask, v1, v2, correct);
        
        //TODO: selectEqual, selectNotEqual...

        static if(isFloatingPoint!T)
        {
            v1 = vector!T(seq);
            v2 = vector!T(staticIota!(2, 2 + n));
            
            auto magSq = (T[] a) => reduce!"a + b * b"(to!T(0), a);
            
            static if(n >= 2)
                test!dot2(v1, v2, repeated!V(1 * 2 + 2 * 3));
            
            static if(n >= 3)
            {
                test!dot3(v1, v2, repeated!V(1 * 2 + 2 * 3 + 3 * 4));
                test!(magnitude3, true)(v1, repeated!V(
                    std.math.sqrt(magSq(v1.ptr_[0 .. 3]))));
                
                test!(magEst3, true)(v1, repeated!V(
                    std.math.sqrt(magSq(v1.ptr_[0 .. 3]))));
                
                test!(magSq3)(v1, repeated!V(magSq(v1.ptr_[0 .. 3])));

                correct = v1;
                correct /= std.math.sqrt(magSq(v1.ptr_[0 .. 3]));
                test!(normalise3, true)(v1, correct);
                test!(normEst3, true)(v1, correct);
            }

            static if(n >= 4)
            {
                test!dot4(v1, v2, repeated!V(1 * 2 + 2 * 3 + 3 * 4 + 4 * 5));
                test!(magnitude4, true)(v1, repeated!V(
                    std.math.sqrt(magSq(v1.ptr_[0 .. 4]))));
                
                test!(magEst4, true)(v1, repeated!V(
                    std.math.sqrt(magSq(v1.ptr_[0 .. 4]))));
                
                test!(magSq4)(v1, repeated!V(magSq(v1.ptr_[0 .. 4])));

                correct = v1;
                correct /= std.math.sqrt(magSq(v1.ptr_[0 .. 4]));
                test!(normalise4, true)(v1, correct);
                test!(normEst4, true)(v1, correct);
            }
            
            static if(n == 4)
                test!cross3(v1, v2, vector!T(
                    v1.ptr_[1] * v2.ptr_[2] - v1.ptr_[2] * v2.ptr_[1],
                    v1.ptr_[2] * v2.ptr_[0] - v1.ptr_[0] * v2.ptr_[2],
                    v1.ptr_[0] * v2.ptr_[1] - v1.ptr_[1] * v2.ptr_[0], 
                    0));
        }

        alias staticIota!(1, vecBytes + 1) bseq;   
        v = cast(V) vector!ubyte(bseq);

        correct = cast(V) vector!ubyte(bseq[1 .. $], 0);
        test!(shiftBytesLeftImmediate, false, 1, SIMDVer.SSE42)(v, correct); 
        
        correct = cast(V) vector!ubyte(bseq[1 .. $], 1);
        test!(rotateBytesLeftImmediate, false, 1, SIMDVer.SSE42)(v, correct);

        correct = cast(V) vector!ubyte(0, bseq[0 .. $ - 1]);
        test!(shiftBytesRightImmediate, false, 1, SIMDVer.SSE42)(v, correct); 
        
        correct = cast(V) vector!ubyte(vecBytes, bseq[0 .. $ - 1]);
        test!(rotateBytesRightImmediate, false, 1, SIMDVer.SSE42)(v, correct);

        v = vector!T(seq);
        
        correct = vector!T(seq[1 .. $], 0);
        test!(shiftElementsLeft, false, 1, SIMDVer.SSE42)(v, correct); 
        
        correct = vector!T(seq[1 .. $], 1);
        test!(rotateElementsLeft, false, 1, SIMDVer.SSE42)(v, correct);

        correct = vector!T(0, seq[0 .. $ - 1]);
        test!(shiftElementsRight, false, 1, SIMDVer.SSE42)(v, correct); 
        
        correct = vector!T(n, seq[0 .. $ - 1]);
        test!(rotateElementsRight, false, 1, SIMDVer.SSE42)(v, correct);

        {
            alias staticIota!(1, 2 * n + 1) both;
            v1 = vector!T(both[0 .. n]);
            v2 = vector!T(both[n .. $]);
    
            correct = vector!T(both[1 .. n + 1]);
            test!(shiftElementsLeftPair, false, 1, SIMDVer.SSE42)(v1, v2, correct); 
            
            correct = vector!T(both[n - 1 .. $ - 1]);
            test!(shiftElementsRightPair, false, 1, SIMDVer.SSE42)(v2, v1, correct); 
        }

        v1 = vector!T(seq);
        v2 = v1;
        v2.ptr_[0] = 2; 
        
        mask = ubyte.max;
        mask.ptr_[0 .. T.sizeof] = 0;
        test!maskEqual(v1, v2, mask);

        mask = 0;
        mask.ptr_[0 .. T.sizeof] = ubyte.max; 
        test!maskNotEqual(v1, v2, mask);
    
        mask = 0;
        mask.ptr_[0 .. T.sizeof] = ubyte.max;
        test!maskGreater(v2, v1, mask);

        mask = ubyte.max;
        mask.ptr_[0 .. T.sizeof] = 0;
        test!maskGreaterEqual(v1, v2, mask);

        mask = 0;
        mask.ptr_[0 .. T.sizeof] = ubyte.max;
        test!maskLess(v1, v2, mask);

        mask = ubyte.max;
        mask.ptr_[0 .. T.sizeof] = 0;
        test!maskLessEqual(v2, v1, mask);
    
        version(none) // not implemented yet
        {
            v3 = vector!T(staticIota!(2, 2 + n));
            
            test!allEqual(v1, v1, true);
            test!allEqual(v1, v2, false);
            test!allNotEqual(v1, v3, true);
            test!allNotEqual(v1, v2, false);
        }
         
    }
    
    {
        auto mask = vector!ubyte(
            10, 4, 0, 1, 2, 11, 3, 6, 7, 8, 2, 9, 12, 6, 14, 15);
        auto v = vector!ubyte(staticIota!(0, 16));
        test!permute(v, mask, mask);
    }

    //TODO: toFloat, toDouble, toInt for cases when the argument and the 
    // result do not have the same number of elements.
}
