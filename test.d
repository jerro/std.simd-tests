import std.stdio, std.traits, std.math;
import std.simd;
import core.simd;

template staticIota(size_t n, r...)
{
    static if(n == 0)
        alias r staticIota;
    else 
        alias staticIota!(n - 1, n - 1, r) staticIota;
}

auto vector(A...)(A a)
{
    alias Vector!(A[0][A.length]) V;
    V r;
    foreach(i,_; A)
        (cast(A[0]*) &r)[i] = a[i];
    
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
        float[2] a;
        storeUnaligned(v, &a[1]);
        assert(equal(a[1], 1f));
    });
    
    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        auto x = vector(5f, 5, 5, 5);
        v = setX(v, x);
        assert(equal(v, vector(5f, 2, 3, 4)));
    });
    
    test!((_)
    {
    });
}
