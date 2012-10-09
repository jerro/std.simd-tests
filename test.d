import std.stdio, std.traits, std.math, std.conv, std.typetuple, std.string;
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

auto vectorT(T, A...)(A a)
{
    alias Vector!(T[A.length]) V;
    V r;
    foreach(i,_; A)
        (cast(T*) &r)[i] = cast(T) a[i];
    
    return r;
}

auto vector(A...)(A a) { return vectorT!(A[0])(a); }

auto eq(T)(T a, T b) if(isIntegral!T)
{
    return a == b;
}

auto eq(T)(T a, T b) if(isFloatingPoint!T)
{
    return feqrel(a, b) + 2 >= T.mant_dig;
}

auto eq(V)(V a, V b) if(!isFloatingPoint!V && !isIntegral!V)
{
    alias typeof(a.array[0]) T;

    foreach(i; staticIota!(0, V.init.length))
        if(!eq(a.array[i], b.array[i]))
            return false;

    return true;
}

void test(alias f)(int line = __LINE__, string message = "")
{
    static if(is(typeof(f(0))))
        f(0);
    else
        writefln("Test on line %s failed to compile.", line); 
}

template group(a...){ alias a members; }

template ungroup(g...) if(g.length == 1)
{
    alias g[0] g0;

    static if(is(typeof(g0.members)))
        alias g0.members ungroup;
    else
        alias g0 ungroup;
}

template testElementWise_(
    alias finit, alias vecOps, alias fresult, int vecBytes = 16, Scal...)
{
    static if(Scal.length == 0)
        alias tt!(float, double, byte, ubyte, short, ushort, int, uint, long, ulong)
            Scalars;
    else 
        alias Scal Scalars; 

    template initGroup(int n)
    {
        template initGroup(alias i)
        {
            alias group!(staticMapF!(ungroup!(finit)[i], staticIota!(0, n))) 
                initGroup;
        }
    }

    void testElementWise()
    {
        foreach(T; Scalars)
        {
            enum n = vecBytes / T.sizeof;
            enum nParams = ungroup!(finit).length;
            alias staticMap!(initGroup!n, staticIota!(0, nParams)) initTuple;

            alias Vector!(T[n]) V;
            staticRepeat!(nParams, V) params;
            foreach(i, _; params)
                params[i] = vectorT!T(ungroup!(initTuple[i])); 
            
            foreach(i, _; params)
                writeln(params[i].array);
 
            writeln();

            foreach(i, _; vecOps)
            {
            }
        }
    } 
}

void testElementWise(
    alias finitGroup, alias opsGroup, int vecBytes = 16, Scal...)()
{
    static if(Scal.length == 0)
        alias tt!(float, double, byte, ubyte, short, ushort, int, uint, long, ulong)
            Scalars;
    else 
        alias Scal Scalars; 
   
    alias ungroup!finitGroup finit;
    alias ungroup!opsGroup ops;
 
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
            
            static if(is(typeof(ops[i](params))))
                assert(eq(ops[i](params), correct), format(
                    "Function %s returned an incorrect result when called with parameters of type %s",
                     ops[i].stringof, V.stringof ));
            else
                pragma(msg, 
                    "Failed to compile function " ~ ops[i].stringof ~ 
                    " called with parameters of type " ~ V.stringof); 
        }
    }
}

void main()
{
    testElementWise!(
        group!(a => 3 * a, a => a + 5), 
        group!(
            add, (a, b) => a + b,  
            mul, (a, b) => a * b))();

    enum vecBytes = 16;
    
    test!((_)
    {
        float a = 1;
        auto v = loadScalar!float4(&a);
        assert(eq(v, vector(1f, 0, 0, 0))); 
    });

    test!((_)
    {
        float[5] a = [0, 1, 2, 3, 4];
        auto v = loadUnaligned!float4(&a[1]);
        assert(eq(v, vector(1f, 2, 3, 4)));
    });
    
    test!((_)
    {
        auto v = vector(1f, 2, 3, 4);
        assert(eq(getScalar(v), 1f));
    });
    
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
        foreach(T; tt!(int4, float4))
        {
            auto v1 = cast(T) vector(1, 2, 3, 4);
            auto v2 = cast(T) vector(1 | 8, 2 | 8, 3 | 8, 4 | 8);
            auto v3 = cast(T) vector(8, 8, 8, 8);
            assert(eq(xor(v1, v2), v3));
        }
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

    test!((_)
    {
        auto v1 = vector(1f, 2, 3, 4);
        auto v2 = vector(2f, 2, 2, 2);
        assert(eq(cast(int4) maskGreater(v1, v2), vector(0, 0, -1, -1)));
    });

    test!((_)
    {
        auto f32 = vector(1f, -2, 3, -4);
        assert(eq(abs(f32), vector(1f, 2, 3, 4))); 
        
        auto i32 = vector(1, -2, 3, -4);
        assert(eq(abs(i32), vector(1, 2, 3, 4))); 
        
        auto f64 = vector(1.0, -2);
        assert(eq(abs(f64), vector(1.0, 2)));
    });

    test!((_)
    {
        foreach(T; tt!(byte, ubyte, short, ushort))
        foreach(saturate; tt!(true, false))
        foreach(addition; tt!(true, false))
        static if(std.traits.isSigned!T || addition)
        {
            enum s = T.max;
            enum n = vecBytes / T.sizeof;        
            auto v1 = vectorT!T(staticIota!(1, 1 + n));
            auto r = vectorT!T(staticIota!(2, 2 + 2 * n, 2));

            static if(saturate)
            { 
                typeof(v1) y = cast(T)(11 + s / 2);
                v1 = setY(v1, y);
                y = cast(T)s; 
                r = setY(r, y);
            }

            static if(addition)
            {
                alias staticSelect!(saturate, addSaturate, add) op;
                auto v2 = v1;
            }
            else
            {
                alias staticSelect!(saturate, subSaturate, sub) op;
                auto v2 = -v1;
            }

            assert(eq(op(v1, v2), r));
        }
    });

    test!((_)
    {
        auto v1 = vector(1f, 2, 3, 4);
        auto v2 = vector(2f, 3, 4, 5);
        auto v3 = vector(3f, 4, 5, 6);

        assert(eq(mul(v1, v2), vector(2f, 6, 12, 20)));
        assert(eq(madd(v1, v2, v3), vector(5f, 10, 17, 26)));
        assert(eq(msub(v1, v2, v3), vector(-1f, 2, 7, 14)));
        assert(eq(nmadd(v1, v2, v3), vector(1f, -2, -7, -14)));
        assert(eq(nmsub(v1, v2, v3), vector(-5f, -10, -17, -26)));
    });
    
    // TODO: selectGreater

    test!((_)
    {
        foreach(T; tt!(float, double, int, uint, short, ushort, byte, ubyte))
        {
            enum n = vecBytes / T.sizeof;
            auto v1 = vectorT!T(staticIota!(1, n + 1));
            auto v2 = vectorT!T(staticRepeat!(n, 2));
            auto v3 = vectorT!T(staticRepeat!(n, 3));
            auto vmin = vectorT!T(1, staticRepeat!(n - 1, 2));
            auto vmax = vectorT!T(2, staticIota!(2, n + 1));
            auto vclamp = vectorT!T(2, 2, staticRepeat!(n - 2, 3));
            assert(eq(min!(SIMDVer.SSE41)(v1, v2), vmin));
            assert(eq(max!(SIMDVer.SSE41)(v1, v2), vmax));
            assert(eq(clamp!(SIMDVer.SSE41)(v2, v1, v3), vclamp));
        }
    });

    test!((_)
    {
        foreach(T; tt!(float, double))
        {
            enum n = vecBytes / T.sizeof;
            auto v1 = vectorT!T(staticRepeat!(n, 1));
            auto v2 = vectorT!T(staticRepeat!(n, 1 + n));
            auto dt = 1.0 / n;
            auto t = vectorT!T(staticIota!(0, n));
            t /= n;
            auto r = vectorT!T(staticIota!(1, n + 1));
            assert(eq(lerp(v1, v2, t), r));
        }
    });

    test!((_)
    {
        foreach(T; tt!(float, double))
        {
            enum n = vecBytes / T.sizeof;
            alias staticMapF!(a => 0.55 + 0.9 * a, staticIota!(0, n)) vtuple;
            auto v = vectorT!T(vtuple);
            auto rfloor = vectorT!T(staticMapF!(a => cast(int) a, vtuple)); 
            auto rceil = vectorT!T(staticMapF!(a => cast(int)(a + 1), vtuple)); 
            auto rround = vectorT!T(staticMapF!(a => cast(int)(a + 0.5), vtuple));
            assert(eq(floor!(SIMDVer.SSE41)(v), rfloor));
            assert(eq(trunc!(SIMDVer.SSE41)(v), rfloor));
            assert(eq(ceil!(SIMDVer.SSE41)(v), rceil));
            assert(eq(round!(SIMDVer.SSE41)(v), rround));
        }
    }); 
}
