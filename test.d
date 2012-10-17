import std.stdio, std.traits, std.math, std.conv, std.typetuple, std.string,
    std.bigint;
import std.simd;
import core.simd;

alias TypeTuple tt;

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

void print(T)(T a)
{
    static if(is(typeof(a.array)))
        writeln(a.array);
    else
        writeln(a);
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
            assert(eq!approx(op!tParams(params), correct), format(
                "Function %s, using instructions set %s,"
                " returned an incorrect result"
                " when called with parameters of type %s",
                op.stringof, tParams[$-1].stringof, typeof(params).stringof));
        }
        else
            pragma(msg, 
                "Failed to compile: " ~ op.stringof ~
                "   parameters: " ~ typeof(params).stringof ~ 
                "   ver:  " ~ tParams[$ - 1].stringof); 
    }
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

auto testStoreScalar(SIMDVer ver, V)(V vec)
{
    BaseType!V s;
    storeScalar!ver(vec, &s);
    return s;
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
 
        T s = 1;
        V correct = 0;
        correct.array[0] = 1;
        test!(loadScalar, false, V, SIMDVer.SSE42)(&s, correct);
       
        T[n + 1] a = [0, seq];
        test!(loadUnaligned, false, V, SIMDVer.SSE42)(&a[1], vector!T(seq));

        testStoreScalar!(SIMDVer.SSE42,)(vector!T(seq));
        test!testStoreScalar(vector!T(seq), to!T(1));
    }

    /* 
    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        float a;
        storeScalar(v, &a);
        assert(eq(a, 1f));
    });

    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        float[5] a;
        storeUnaligned(v, &a[1]);
        assert(eq(a[1], 1f));
    });

    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        assert(eq(getScalar(v), 1f));
    });
    
    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        auto x = vector(5f, 5, 5, 5);
        assert(eq(setX(v, x), vector(5f, 2, 3, 4)));
    });
    
    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        auto y = vector(5f, 5, 5, 5);
        assert(eq(setY(v, y), vector(1f, 5, 3, 4)));
    });

    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        auto z = vector(5f, 5, 5, 5);
        assert(eq(setZ(v, z), vector(1f, 2, 5, 4)));
    });

    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        auto w = vector(5f, 5, 5, 5);
        assert(eq(setW(v, w), vector(1f, 2, 3, 5)));
    });
    
    test!((_)
    {
        foreach(i; staticIota!(0, 256))
        {
            alias tt!(i & 3, (i >> 2) & 3, (i >> 4) & 3, (i >> 6) & 3) e;
            auto v = vector(0f, 1, 2, 3);
            v = swizzle!(concatenate!e)(v);
            assert(eq(v, vector(staticMap!(toFloatTemplate, e)))); 
        } 
    });

    test!((_)
    {
        auto v = vector(staticMap!(toUByte, staticIota!(0, vecBytes)));
        auto mask = vector(
            cast(ubyte) 10, 4, 0, 1, 2, 11, 3, 6, 7, 8, 2, 9, 12, 6, 14, 15); 

        v = permute!(SIMDVer.SSE3)(v, mask);
        assert(eq(v, mask));
    });

    test!((_)
    {
        auto v1 = vector(0f, 1, 2, 3);
        auto v2 = vector(4f, 5, 6, 7);
        assert(eq(interleaveLow(v1, v2), vector(0f, 4, 1, 5))); 
        assert(eq(interleaveHigh(v1, v2), vector(2f, 6, 3, 7))); 
    });
    
    test!((_)
    {
        auto v = vector(cast(short) 0, 1, 2, 3, 4, 5, 6, 7);
        assert(eq(unpackLow(v), vector(0, 1, 2, 3))); 
        assert(eq(unpackHigh(v), vector(4, 5, 6, 7))); 
    });
    
    test!((_)
    {
        auto v1 = vector(0, 1, 2, 3);
        auto v2 = vector(4, 5, 6, 7);
        assert(eq(pack(v1, v2), vector(cast(short) 0, 1, 2, 3, 4, 5, 6, 7)));
    });
    
    test!((_)
    {
        auto s = short.max;
        auto v1 = vector(0, 1 + s, 2, 3);
        auto v2 = vector(4, 5, 6 + s, 7);
        assert(eq(packSaturate(v1, v2), 
            vector(cast(short) 0, s, 2, 3, 4, 5, s, 7)));
    });

    test!((_)
    {
        assert(eq(toInt(vector(0f, 1, 2, 3)), vector(0, 1, 2, 3))); 
        assert(eq(toInt(vector(1.0 , 2)), vector(1, 2, 0, 0))); 
    });

    test!((_)
    {
        assert(eq(toFloat(vector(0, 1, 2, 3)), vector(0f, 1, 2, 3))); 
        assert(eq(toFloat(vector(1.0 , 2)), vector(1f, 2, 0, 0))); 
    });

    test!((_)
    {
        assert(eq(toDouble(vector(1, 2, 3, 4)), vector(1.0, 2))); 
        assert(eq(toDouble(vector(1f, 2, 3, 4)), vector(1.0, 2))); 
    });

    test!((_)
    {
        auto v = vector(1 << 3, 2 << 3, 3 << 3, 4 << 3);
        assert(eq(shiftRightImmediate!3(v), vector(1, 2, 3, 4))); 
    });
    
    test!((_)
    {
        foreach(T; tt!(int, float))
        {
            auto v1 = vectorT!T(1, 2, 3, 4);
            auto v2 = vectorT!T(5, 6, 7, 8);
            auto mask = cast(void16) vector(0, -1, 0, -1);
            assert(eq(std.simd.select(mask, v1, v2), vectorT!T(1, 6, 3, 8)));
        } 
    });
    */
}
