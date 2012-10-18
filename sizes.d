import std.stdio, std.process, std.string, std.range, std.algorithm, std.conv;

void main()
{
    bool[string] table;
    foreach(l; std.algorithm.splitter(shell("./test -l")))
        if(!l.canFind("test"))
            table[l] = true;

    struct S{ string name; size_t size; }
    S[] arr;
    
    foreach(l; std.algorithm.splitter(shell("nm -St decimal test"), "\n"))
    {
        auto s = split(l);
        if(s.length != 4)
            continue;

        auto name = s.back;
        if(name in table)
            arr ~= S(name, to!size_t(s[1]));
    }

    sort!"a.size < b.size"(arr);
    foreach(e; arr)
        writefln("%s %s", e.tupleof);
}
