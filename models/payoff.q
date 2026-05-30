/ payoff.q — terminal payoff calculation

.payoff.calculateVanillaPayoff:{[optType;spotGrid;strike]
    $[optType~`call; 0f|spotGrid-strike;
      optType~`put;  0f|strike-spotGrid;
      '"Unsupported optionType for payoff: ",string optType]
 };

.payoff.calculateIntrinsicValue:{[trade;spotGrid]
    .payoff.calculateVanillaPayoff[trade`optionType;spotGrid;trade`strike]
 };
