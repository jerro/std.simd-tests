import std.stdio, std.traits, std.math, std.conv, std.typetuple, std.string,
    std.bigint, std.random, std.algorithm, std.range;
import std.simd;
import core.simd;

alias TypeTuple tt;

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

template staticMapF(alias f, a...)
{
    static if(a.length == 0)
        alias tt!() staticMapF;
    else
        alias tt!(f(a[0]), staticMapF!(f, a[1 .. $])) staticMapF;
}

template concatenate(a...)
{
    static if(a.length == 0)
        enum concatenate = "";
    else 
        enum concatenate = to!string(a[0]) ~ concatenate!(a[1 .. $]);
}

template staticSelect(bool first, a...) if (a.length == 2)
{
    static if(first)
        alias a[0] staticSelect;
    else
        alias a[1] staticSelect; 
}

template toFloatTemplate(alias a)
{
    enum toFloatTemplate = cast(float) a;
}

template toUByte(alias a)
{
    enum toUByte = cast(ubyte) a;
}

auto vector(T, A...)(A a)
{
    alias Vector!(T[A.length]) V;
    V r;
    foreach(i,_; A)
        (cast(T*) &r)[i] = cast(T) a[i];
    
    return r;
}

auto eq(bool approx = false, T)(T a, T b) if(isIntegral!T)
{
    return a == b;
}

auto eq(bool approx = false, T)(T a, T b) if(isFloatingPoint!T)
{
    if(a.isnan && b.isnan)
        return true;

    return feqrel(a, b) + 3 >= (approx ? T.mant_dig / 2 : T.mant_dig);
}

auto eq(bool approx = false, V)(V a, V b) if(!isFloatingPoint!V && !isIntegral!V)
{
    alias typeof(a.array[0]) T;

    foreach(i; staticIota!(0, V.init.length))
        if(!eq!approx(a.array[i], b.array[i]))
            return false;

    return true;
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

template VecTypes(int vecSize)
{
    template VecTypes(A...)
    {
        static if(A.length == 0)
            alias tt!() VecTypes;
        else
            alias tt!(Vector!(A[0][vecSize / A[0].sizeof]), VecTypes!(A[1 .. $])) 
                VecTypes;
    }
}

template BaseType(V)
{
    alias typeof(V.array[0]) BaseType;
}

template nElements(V)
{
    enum nElements = V.sizeof / BaseType!(V).sizeof;
}

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
                
                stderr.writeln("Result: ");
                stderr.print(r);
                
                stderr.writeln("Correct result: ");
                stderr.print(correct);
 
                assert(false, format(
                    "Function %s, using instructions set %s,"
                    " returned an incorrect result"
                    " when called with parameters of type %s",
                    op.stringof, tParams[$-1].stringof, 
                    typeof(params).stringof));
            }           
        }
        else
            pragma(msg, 
                "Failed to compile: " ~ op.stringof ~
                "   parameters: " ~ typeof(params).stringof ~ 
                "   ver:  " ~ tParams[$ - 1].stringof); 
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

auto simpleSwizzle(V)(string ind, V v)
{
    V r;
    foreach(i, char c; ind)
        r.array[i] = v.array["0123456789ABCDEF".countUntil(c)]; 
    
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
    r.array[] = mem[1 .. $];
    return r;
}

void main()
{
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
            std.simd.sqrt, std.math.sqrt),
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
            addSaturate, opSaturate!"+",
            subSaturate, opSaturate!"-",
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

    alias VecTypes!16 Vec;

    foreach(V; Vec!(floatingPoint, integral))
    {
        enum n = nElements!V; 
        alias BaseType!V T;       
        alias staticIota!(1, n + 1) seq;
        
        Vector!(ubyte[V.sizeof]) mask;
        V v, v1, v2; 
        V ones = 1;
        V zeros = 0;
        V correct;
        V vseq = vector!T(seq);

        T s = 1;
        correct = 0;
        correct.array[0] = 1;
        test!(loadScalar, false, V, SIMDVer.SSE42)(&s, correct);
       
        T[n + 1] a = [0, seq];
        test!(loadUnaligned, false, V, SIMDVer.SSE42)(&a[1], vector!T(seq));

        test!getScalar(vector!T(seq), to!T(1));
        test!testStoreScalar(vector!T(seq), to!T(1));
        test!testStoreUnaligned(vector!T(seq), vector!T(seq));

        correct = 0;
        correct.array[0] = 1;
        test!setX(zeros, ones, correct);
        test!getX(correct, ones);

        static if(n >= 2)
        {
            correct = 0;
            correct.array[1] = 1;
            test!setY(zeros, ones, correct);
            test!getY(correct, ones);
        }
        
        static if(n >= 3)
        {
            correct = 0;
            correct.array[2] = 1;
            test!setZ(zeros, ones, correct);
            test!getZ(correct, ones);
        }

        static if(n >= 4)
        {
            correct = 0;
            correct.array[3] = 1;
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
            v.array[i] = 1;
            test!(swizzle, false, swiz, SIMDVer.SSE42)(v, ones); 
        }

        {
            alias staticIota!(0, 2 * n) both;
            v1 = vector!T(both[0 .. n]);
            v2 = vector!T(both[n .. 2 * n]);
            auto result = new T[2 * n];
            foreach(i; both)
               result[i] =  n * (i & 1) + (i >> 1);
            correct.array[] = result[0 .. n];
            test!interleaveLow(v1, v2, correct);
            correct.array[] = result[n .. $];
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
            v2.array[0] = 1;

            test!shiftLeft(v1, v2, vector!T(staticIota!(0, 4 * n, 4))); 
            test!(shiftLeftImmediate, false, 1, SIMDVer.SSE42)(
                v1, vector!T(staticIota!(0, 4 * n, 4))); 
            
            test!shiftRight(v1, v2, vector!T(staticIota!(0, n))); 
            test!(shiftRightImmediate, false, 1, SIMDVer.SSE42)(
                v1, vector!T(staticIota!(0, n))); 
        } 

        mask = 0;
        foreach(i; 0 .. T.sizeof)
            mask.array[i] = ubyte.max;
        v1 = vector!T(staticIota!(50, 50 + n));
        v2 = vector!T(seq);
        correct = v2;
        correct.array[0] = v1.array[0];
        test!(std.simd.select)(mask, v1, v2, correct);
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
