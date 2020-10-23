module Analyzer
using Plots
using Printf
import PyPlot

export analyze

pyplot()
default(show = true)
PyPlot.clf()
println(PyPlot.backend)

function analyze()
    println("Bellaaa!")
    s = "./tags.txt"
    a = readlines(s)
    for y in a
        filter(x -> !isspace(x), y)
    end
    machine_time = 81e-12
    i=0
    b = Array{Int, 2}(undef, 2, length(a))

    b[1, :] = [parse(Int, split(x, ";")[1]) for x in a]
    b[2, :] = [parse(Int, split(x, ";")[2]) for x in a]


    tags = Array{Int, 2}(undef, 3, length(b))
    fill!(tags, 0)
    println(typeof(tags))
    k = Array{Int, 1}(undef, 3) # k[i] will be the total count of trigger events on channel i
    fill!(k, 1)
    i=0
    cnt = 0
    for i = 1:length(a)
        # if (i>7)
        #     println("condition: ", b[1, i] ," < ", tags[ b[2, i]-1, k[b[2, i]-1] - 1 ], "  ?")
        # end
        if (i<8 || b[1, i] - 3900 > tags[ b[2, i]-1, k[b[2, i]-1] - 1 ])
            tags[ b[2, i]-1, k[b[2, i]-1] ] = b[1, i]
            k[b[2, i] - 1] += 1
        else
            println("Afterpulse on CH-", b[2, i] - 1)
        end
    end
    
    println("First 100 triggers on channel 2: ", tags[1, 1:10])

    println("Plotting arrival differences to nearest trigger event")
    g1 = tags[3, 1]
    g2 = tags[3, 2]
    diff1 = Array{Int, 1}(undef, k[1])
    diff2 = Array{Int, 1}(undef, k[2])
    fill!(diff1, 0)
    fill!(diff2, 0)
    n = 1
    for i = 3:k[3]
        while (tags[1, n]<g2 && n<k[1])
            diff1[n] = minimum([abs(tags[1, n] - g1), abs(g2 - tags[1, n])])
            n += 1
        end
        g1 = g2
        g2 = tags[3, i]
    end
    println("maximum difference CH1: ", maximum(diff1)*81e-12)
    max_diff1 = maximum(diff1)

    sort!(diff1)
    diff_advance = 0
    diff_bin_size = 1000000
    grouped = Array{Int, 1}(undef, Int(ceil((max_diff1/diff_bin_size))) + 1)
    fill!(grouped, 0)
    i = 1
    j = 1
    while (i<k[1])
        count = 0
        while (diff1[i] <= diff_advance && i<k[1])
            count += 1
            i += 1
        end
        diff_advance += diff_bin_size
        grouped[j] = count
        j += 1
    end

    println("maximum diff count: ", maximum(grouped))
    hist = Array{Int}(undef, maximum(grouped) + 1)
    fill!(hist, 0)
    i = 1
    while (i<=length(grouped))
        hist[grouped[i] + 1] += 1
        i += 1
    end
    println("Plotting...")
    first_n = length(hist)
    fig = PyPlot.figure()
    fig = Plots.bar(1:first_n, hist[1:first_n], show=true)
    display(fig)
    return 0
end

function poisson_pdf(x, mu)
    return (mu)^x/factorial(big(x)) * exp(-mu)
end

function bose_einst_pdf(x, mu)
    return 1/(mu + 1) * (mu/(mu+1))^x
end

function poisson_moments(mu)
    return [mu, mu, mu, mu]
end

function bose_ein_moments(mu)
    return [mu, mu + mu^2, mu + 3*mu^2 + 2*mu^3, mu + 10*mu^2 + 18*mu^3 + 9*mu^4]
end

function main()
    inter = 20e-6
    analyze()
end
end