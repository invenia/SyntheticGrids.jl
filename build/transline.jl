@auto_hash_equals struct TransLine
    connecting::Tuple{Bus,Bus}
    impedance::Real # Ohms
    capacity::Real # A

    function TransLine(
        connecting::Tuple{Bus,Bus};
        impedfunc=linear_imped,
        capfunc=volt_cap
    )
        return new(
            connecting,
            impedfunc(connecting[1], connecting[2]),
            capfunc(connecting[1], connecting[2])
        )
    end
end

function TransLine(a::Bus, b::Bus; impedfunc=linear_imped,  capfunc=volt_cap)
    return TransLine((a, b), impedfunc=impedfunc, capfunc=capfunc)
end

function show(io::IO, tl::TransLine)
    println(io, "TransLine(")
    println(io, "\tconnecting: $(tl.connecting),")
    println(io, "\timpedance=$(tl.impedance),")
    println(io, "\tcapacity=$(tl.capacity)")
    print(io, ")")
end

"""
    linear_imped(a::Bus, b::Bus)

Return an impedance value for a transmission line connecting buses `a` and `b`, using a
simple linear approximation. This value is loosely based on the two references below.
Values based on resistances for dc currents at 50C and 350kV lines.

REF 1: Electric Power Research Institute (EPRI).
Transmission Line Reference Book, 345 kV and Above.
Palo Alto, CA: The Institute, 1975.
REF 2: J. Glover, M. Sarma and T. Overbye, Power System Analysis and Design.
Stamford, CT: Cengage Learning, 2012.
"""
function linear_imped(a::Bus, b::Bus)
    respkm = 0.045 # Ohms/km
    d = haversine(a, b) # this is using the
    # straight line distance between nodes. Considering the lines arc, real length
    # is slightly longer. However, this average value is already approximated
    # enough.
    return d * respkm
end

"""
    rand_imped(a::Bus, b::Bus)

Return an impedance value for a transmission line connecting buses `a` and `b`, using a
simple linear approximation with a randomly chosen coefficient. The range we adopted is
based on the values for ACSR cables from the ref: J. Glover, M. Sarma and T. Overbye, Power
System Analysis and Design. Stanford, CT: Cengage Learning, 2012.
"""
function rand_imped(a::Bus, b::Bus)
    respkm = 0.025:0.001:0.250
    d = haversine(a, b)
    return d * rand(respkm)
end

"""
    volt_cap(a::Bus, b::Bus)

Return a current carrying capacity value for a transmission line connecting buses `a` and
`b` based on their operating voltages. This value is loosely based on the reference:
Electric Power Research Institute (EPRI). Transmission Line Reference Book, 345 kV and
Above. Palo Alto, CA: The Institute, 1975.
"""
function volt_cap(a::Bus, b::Bus)
    cappvolt = 14 # A/kV
    volt1 = a.voltage
    volt2 = a.voltage
    bvolt = [v for v in volt1 if v in volt2]
    svolt = vcat(a.voltage, b.voltage)
    if length(bvolt) > 0
        return bvolt[rand(1:length(bvolt))] * cappvolt
    else
        return svolt[rand(1:length(svolt))] * cappvolt
    end
end

end_points(t::TransLine) = t.connecting
impedance(t::TransLine) = t.impedance
capacity(t::TransLine) = t.capacity
