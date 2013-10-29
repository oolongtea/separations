- heatXchanger.m
> returns cost of heat exchanger 
> inputs: heat load [MW], approach temperatures [C], overall heat coefficient, material factors

- pump.m
> return cost of pump(s)
- inputs: (array of) amount to be pumped [g/s], pressure [torr]

- distillation_column.m
> script that prints conditions and compositions of a heptane-pentance distillation column

/BeerBrewModel/
- beer_brew.m
> script that uses a beer fermentation model that predicts the final alcohol content by alcohol by volume (ABV), time for full fermenation, and size of the heat exchanger needed to mainatin isothermal conditions (using the Monod equations)
>> muone.m, mutwo.m, muthree.m are the Monod equations
