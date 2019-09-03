#include <getopt.h>
#include <iostream>
#include <fstream>
#include <cmath>

// May use std::optional in C++ 17
// https://stackoverflow.com/questions/6822044/checking-if-a-variable-is-initialized

double T,ystart, dt;
std::string outPath = "-";

void PrintHelp()
{
    std::cout <<
        "Parameter settings:\n"
        "   -T, --maxTime:         Set maximum time\n"
        "   -y, --initialState:    Set initial condition\n"
        "   -d, --timeStep:        Set time step delta t\n"
        "   -o, --outfile:         File to write solution\n"
        "   -h, --help:            Show help\n";
    exit(1);
}

// https://gist.github.com/ashwin/d88184923c7161d368a9
void ProcessArgs(int argc, char** argv)
{
    const char* const short_opts = "T:y:d:o:h";
    const option long_opts[] = {
        {"maxTime", required_argument, nullptr, 'T'},
        {"initialState", no_argument, nullptr, 'y'},
        {"timeStep", required_argument, nullptr, 'd'},
        {"outfile", required_argument, nullptr, 'o'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, no_argument, nullptr, 0}
    };

    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt)
        {
        case 'T':
            T = std::stoi(optarg);
            std::cout << "Num set to: " << T << std::endl;
            break;

        case 'y':
            ystart = std::stof(optarg);
            std::cout << "Beep is set to true\n";
            break;

        case 'd':
            dt = std::stof(optarg);
            std::cout << "Sigma set to: " << dt << std::endl;
            break;

        case 'o':
            outPath = std::string(optarg);
            std::cout << "Write file set to: " << outPath << std::endl;
            break;

        case 'h': // -h or --help
        case '?': // Unrecognized option
        default:
            PrintHelp();
            break;
        }
    }
}

double dydt(double y, double t) {
    return sqrt(y);
}

double euler(double y0, double dt, double T, std::ostream& out = std::cout) {
    
    double y = y0;
    double t = 0;
    
    // Write initial state
    out << t << "\t" << y0 << std::endl;
    while (t < T) {
        // Reduce dt such that T is exactly reached. 
        if (t + dt > T) {
            dt = T - t;
        }
        // Evalute for t + dt
        y = y + dt*dydt(y, t);
        t += dt;
        out << t << "\t" << y << std::endl;
    }
    return y;
}

int main(int argc, char* argv[])
{
    ProcessArgs(argc, argv);

    // https://stackoverflow.com/questions/428630/assigning-cout-to-a-variable-name
    std::ofstream outFile;
    if (outPath != "-") {
        outFile.open(outPath, std::ios::out);
    } 
    std::ostream& out = (outPath != "-" ? outFile : std::cout);
    
    euler(ystart, dt, T, out);

    return 0;
}

