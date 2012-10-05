import std.stdio, std.traits, std.math, std.conv, std.typetuple;
import std.simd;
import core.simd;

template staticIota(size_t n, r...)
{
    static if(n == 0)
        alias r staticIota;
    else 
        alias staticIota!(n - 1, n - 1, r) staticIota;
}

template concatenate(a...)
{
    static if(a.length == 0)
        enum concatenate = "";
    else 
        enum concatenate = to!string(a[0]) ~ concatenate!(a[1 .. $]);
}

template toFloatTemplate(alias a)
{
    enum toFloatTemplate = cast(float) a;
}

template toUByte(alias a)
{
    enum toUByte = cast(ubyte) a;
}

auto vector(A...)(A a)
{
    alias Vector!(A[0][A.length]) V;
    V r;
    foreach(i,_; A)
        (cast(A[0]*) &r)[i] = cast(A[0]) a[i];
    
    return r;
}

auto equal(T)(T a, T b) if(isIntegral!T)
{
    return a == b;
}

auto equal(T)(T a, T b) if(isFloatingPoint!T)
{
    return feqrel(a, b) + 2 >= T.mant_dig;
}

auto equal(V)(V a, V b) if(!isFloatingPoint!V && !isIntegral!V)
{
    alias typeof(a.array[0]) T;

    foreach(i; staticIota!(V.init.length))
        if(!equal(a.array[i], b.array[i]))
            return false;

    return true;
}

void test(alias f)(int line = __LINE__)
{
    static if(is(typeof(f(0))))
        f(0);
    else
        writefln("Test on line %s failed to compile.", line); 
}

template ArrayType(T : T[]) { alias T ArrayType; }

void main()
{
    test!((_)
    {
        float a = 1;
        auto v = loadScalar!float4(&a);
        assert(equal(v, vector(1f, 0, 0, 0))); 
    });

    test!((_)
    {
        float[5] a = [0, 1, 2, 3, 4];
        auto v = loadUnaligned!float4(&a[1]);
        assert(equal(v, vector(1f, 2, 3, 4)));
    });
    
    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        assert(equal(getScalar(v), 1f));
    });
    
    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        float a;
        storeScalar(v, &a);
        assert(equal(a, 1f));
    });

    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        float[5] a;
        storeUnaligned(v, &a[1]);
        assert(equal(a[1], 1f));
    });

    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        auto x = vector(5f, 5, 5, 5);
        assert(equal(setX(v, x), vector(5f, 2, 3, 4)));
    });
    
    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        auto y = vector(5f, 5, 5, 5);
        assert(equal(setY(v, y), vector(1f, 5, 3, 4)));
    });

    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        auto z = vector(5f, 5, 5, 5);
        assert(equal(setZ(v, z), vector(1f, 2, 5, 4)));
    });

    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        auto w = vector(5f, 5, 5, 5);
        assert(equal(setW(v, w), vector(1f, 2, 3, 5)));
    });
    
    test!((_)
    {
        foreach(i; staticIota!256)
        {
            alias TypeTuple!(i & 3, (i >> 2) & 3, (i >> 4) & 3, (i >> 6) & 3) e;
            auto v = vector(0f, 1, 2, 3);
            v = swizzle!(concatenate!e)(v);
            assert(equal(v, vector(staticMap!(toFloatTemplate, e)))); 
        } 
    });

    test!((_)
    {
        auto v = vector(staticMap!(toUByte, staticIota!16));
        auto mask = vector(
            cast(ubyte) 10, 4, 0, 1, 2, 11, 3, 6, 7, 8, 2, 9, 12, 6, 14, 15); 

        v = permute!(SIMDVer.SSE3)(v, mask);
        assert(equal(v, mask));
    });

    test!((_)
    {
        auto v1 = vector(0f, 1, 2, 3);
        auto v2 = vector(4f, 5, 6, 7);
        assert(equal(interleaveLow(v1, v2), vector(0f, 4, 1, 5))); 
        assert(equal(interleaveHigh(v1, v2), vector(2f, 6, 3, 7))); 
    });
    
    test!((_)
    {
        auto v = vector(cast(short) 0, 1, 2, 3, 4, 5, 6, 7);
        assert(equal(unpackLow(v), vector(0, 1, 2, 3))); 
        assert(equal(unpackHigh(v), vector(4, 5, 6, 7))); 
    });
    
    test!((_)
    {
        auto v1 = vector(0, 1, 2, 3);
        auto v2 = vector(4, 5, 6, 7);
        assert(equal(pack(v1, v2), vector(cast(short) 0, 1, 2, 3, 4, 5, 6, 7)));
    });
    
    test!((_)
    {
        auto s = short.max;
        auto v1 = vector(0, 1 + s, 2, 3);
        auto v2 = vector(4, 5, 6 + s, 7);
        assert(equal(packSaturate(v1, v2), 
            vector(cast(short) 0, s, 2, 3, 4, 5, s, 7)));
    });

    test!((_)
    {
        assert(equal(toInt(vector(0f, 1, 2, 3)), vector(0, 1, 2, 3))); 
        assert(equal(toInt(vector(1.0 , 2)), vector(1, 2, 0, 0))); 
    });

    test!((_)
    {
        assert(equal(toFloat(vector(0, 1, 2, 3)), vector(0f, 1, 2, 3))); 
        assert(equal(toFloat(vector(1.0 , 2)), vector(1f, 2, 0, 0))); 
    });

    test!((_)
    {
        assert(equal(toDouble(vector(1, 2, 3, 4)), vector(1.0, 2))); 
        assert(equal(toDouble(vector(1f, 2, 3, 4)), vector(1.0, 2))); 
    });
}
