#' Function to write a random completion message.
#' @export
#' @return None
print_completion_message = function(){
  # We have a list with ten random completion messages.
  messages = list()
  messages[[1]]= "Cool beans friend! Thanks for using KBoost."
  messages[[2]]= "KBoost has finished KBoosting your data! Thanks for using KBoost amig@!."
  messages[[3]]= "You sure sound smart for using KBoost! Thanks!"
  messages[[4]]= "KBoost has finished running.(We're feeling serious today)"
  messages[[5]]= "Gracias por usar KBoost! (Feeling Spanish today :D)"
  messages[[6]]= "Ihr Netzwerk ist bereit. Danke.(Feeling German today :D)"
  messages[[7]]= "La sua GRN e prontissima! Grazie per utilizare KBoost.(Feeling Italian today :D)"
  messages[[8]]= "La GRN esta a punt! Moltes gracies per utilizar KBoost. (Feeling Catalan today :D)"
  messages[[9]]=  "Fardig! Takk Takk! (Feeling Swedish today :D)"
  messages[[10]]= "Klaar! dankjewel vriend! (Feeling Dutch today :D)"
  messages[[11]] = "Ta ya lista la tu GRN. Muches gracies por usar KBoost! (Feeling Asturian today :D)"
  messages[[12]] = "Go raibh math agat! (Feeling Irish today :D) "
  messages[[13]] = "Oh My God Hun! You're such a star for using KBoost like I can't cope. xoxo. (We're feeling pretty today :D)"
  # Randomly choose a message.
  idx =  sample(1:length(messages),1)
  print(messages[[idx]])
}
