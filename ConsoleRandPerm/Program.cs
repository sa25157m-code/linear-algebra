using System.CommandLine;
using System.CommandLine.Parsing;
using static Statistics.RandomSampling;

internal class Program
{
    private static int Main(string[] args)
    {
        Option<int> sizeOption = new("--size")
        {
            Description = "Size of the random permutation",
            DefaultValueFactory = parseResult => 0,
            Required = false,
            Arity = ArgumentArity.ExactlyOne
        };

        RootCommand rootCommand = new("ConsoleRandPerm");
        rootCommand.Options.Add(sizeOption);

        ParseResult parseResult = rootCommand.Parse(args);
        if (parseResult.Errors.Count > 0)
        {
            foreach (ParseError parseError in parseResult.Errors)
            {
                Console.Error.WriteLine(parseError.Message);
            }
            return 1;
        }

        var s = Math.Max(5, parseResult.GetValue(sizeOption));

        var perm = RandPerm(s);

        System.Console.WriteLine($"Random permutation of size {s}:");
        for (int i = 0; i < s - 1; i++)
        {
            Console.Write($"{perm[i]},");
        }
        Console.WriteLine($"{perm[s - 1]}");

        return 0;
    }
}