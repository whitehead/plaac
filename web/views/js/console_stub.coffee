# REVIEW: thank god this bit is done.  screw missing console problems.  template.
if window.console is undefined

  $div_console = $("<div id='console' style='display:none;'></div>")

  log = (type,msg) ->
    $div_console.append("<div class='console #{type}'>#{type}: #{msg}</div>")

  $ =>
    $('body').append($div_console)
    $('body').dblclick (e) =>
      if e.ctrlKey
        $div_console.toggle()

  window.console = {
    dir:   (o...) => log('dir',o)
    error: (o...) => log('error',o)
    warn:  (o...) => log('warn',o)
    info:  (o...) => log('info',o)
    log:   (o...) => log('log',o)
    debug: (o...) => log('debug',o)
    trace: (o...) => log('trace',o)
  }

window.p = (msg) ->
  console.debug msg

# vim: ft=coffee
