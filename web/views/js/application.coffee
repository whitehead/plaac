class PickList
  constructor: ->

    # ui
    @picklist_id = '#picklist'
    @candidates_id = '#candidates'

    # cookie
    @cookie_name = 'plaac_finder_picklist'

    # initialization
    @load()
    @render()

  # API

  add: (idx) ->
    #console.log "add #{idx}"
    @picklist.push parseInt(idx)
    @save()
    @picklist_ui_add(idx)
    @render()

  remove: (idx) ->
    #console.log "remove #{idx}"
    arr_idx = _.indexOf(@picklist, parseInt(idx))
    #console.log arr_idx
    if arr_idx >= 0
      @picklist.splice(arr_idx,1)
    @save()
    @render()

  clear: ->
    @picklist = []
    @save()
    @render()

  # UI UPDATES

  picklist_ui_add: (idx) ->
    #console.log "picklist_ui_add: (#{idx})"
    #console.log("#{@candidates_id} #chk#{idx}")
    $chk = $("#{@candidates_id} #chk#{idx}")

    # clone & update ids
    $tr = $chk.parent().parent().clone()
    #console.log($chk)
    #console.log($tr)
    $tr.find('input').attr('id',"picklist_chk#{idx}")
    $tr.find('label').attr('for',"picklist_chk#{idx}")

    title = $tr.find('td:eq(1)').text()
    lbl = "<b>"+title.replace(/\|/,"</b><br/><small>") + "</small>"


    #console.log($tr)
    $("#{@picklist_id}").append """
    <li>
      <label class='checkbox' for='picklist_chk#{idx}' title="#{title}">
        <input type='checkbox' id='picklist_chk#{idx}' name='selected[]' value="#{idx}" checked='checked' />
        #{lbl}
      </label>
    </li>
    """

    ## since we add the element above, 
    ## scroll to keep everything in the same place.
    #$(window).scrollTop($(window).scrollTop() + $tr.height())

  render: ->
    # update picklist
    for chk in $(@picklist_id).find('input[type=checkbox]:not(.all_toggle)')
      $chk = $(chk)
      idx = _.indexOf(@picklist,parseInt($chk.val()))
      is_checked = (idx >= 0)
      $chk.attr('checked',is_checked)
      unless is_checked
        $tr = $chk.parent().parent()
        $tr.remove()

    # update form
    for chk in $(@candidates_id).find('input[type=checkbox]:not(.all_toggle)')
      $chk = $(chk)
      idx = _.indexOf(@picklist,parseInt($chk.val()))
      is_checked = (idx >= 0)
      $chk.attr('checked',is_checked)

    if @picklist.length > 0
      $('#visualize_btn').addClass('btn-primary').removeClass('disabled')
    else
      $('#visualize_btn').removeClass('btn-primary').addClass('disabled')

  # COOKIE IO

  load: ->
    @picklist = [] # default
    cookie_picklist = $.cookie(@cookie_name)
    unless cookie_picklist is null || cookie_picklist == ""
      @picklist = eval(cookie_picklist)

  save: ->
    #console.log('picklist='+@picklist) unless console is undefined
    $.cookie(@cookie_name,"[#{@picklist}]", {path: '/'})
    null

bind_picklist_clear = ->
  $('#run_analysis').click ->
    window.picklist.clear()
  $('.clear_picklist').click ->
    window.picklist.clear()

bind_picklist_toggle = ->
  $('input[type=checkbox]:not(.all_toggle)').click ->
    chk = $(this)
    is_checked = chk.attr('checked')=='checked'
    if is_checked
      picklist.add(chk.val())
    else
      picklist.remove(chk.val())

bind_picklist_toggle_all = ->
  #console.debug "toggle_all"
  $('#candidate_list input[type=checkbox].all_toggle').click ->
    $chk = $(this)
    is_checked = $chk.attr('checked')=='checked'

    $('#candidate_list input[type=checkbox]:not(.all_toggle)').each ->
      chk = $(this)
      chk.attr('checked',is_checked)
      if is_checked
        #console.log "add #{chk.val()}"
        picklist.add(chk.val())
      else
        #console.log "remove #{chk.val()}"
        picklist.remove(chk.val())
      true

bind_hover_highlight = ->

  $('#candidate_list label').mouseover ->
    $(this).parent().parent().addClass('hover')
    $("##{$(this).attr('for')}").parent().parent().addClass('hover')

  $('#candidate_list label').mouseout ->
    $(this).parent().parent().removeClass('hover')
    $("##{$(this).attr('for')}").parent().parent().removeClass('hover')

window.clear_fasta_file = ->
  $('#fasta_file').val('')
  $('#clear_file').hide()
  return `void(0)`

$ ->
  window.picklist = new PickList()

  bind_picklist_clear()
  bind_picklist_toggle()
  bind_picklist_toggle_all()

  bind_hover_highlight()

  handle_slider_change = ->
    bg_freq_val = $('#organisms').val()
    alpha_val = $('#alpha').val()
    if alpha_val is '100'
      $('#organisms').prop('disabled',true)
      $('#organism_background').css('color','#ccc')
      $('.background_probabilities_label').css('color','#ccc')
    else
      $('#organisms').prop('disabled',false)
      $('#organism_background').css('color','black')
      $('.background_probabilities_label').css('color','black')

    if bg_freq_val is '' && alpha_val != '100'
      $('#run_analysis').removeClass('btn-primary')
      $('#run_analysis').prop('disabled',true)
      $('#organisms').css('color','red')
    else
      $('#run_analysis').addClass('btn-primary')
      $('#run_analysis').prop('disabled',false)
      $('#organisms').css('color','black')


  $('#alphaSlider').slider({
    min: -1
    max: 101
    value: 100
    slide: (evt, ui) ->
      $('#alpha').val( $('#alphaSlider').slider('values',0) )
      handle_slider_change()
    change:
      handle_slider_change()
  })

  $('#organisms').on 'change', handle_slider_change
  
  handle_file_change = (e) ->
    val= $(e.target).val()
    # show the clear link if we have a value, else hide it.
    $('#clear_file').toggle(val != "")

    handle_slider_change()

  $('input[name=file]').on 'change', handle_file_change
  handle_file_change({target: $('#fasta_file')})

  # bootstrap
  $('.collapse').collapse(toggle: false)

  $('.collapse').on 'show', ->
    $(this).prev('span').find('i').toggleClass('icon-chevron-right').toggleClass('icon-chevron-down')

  $('.collapse').on 'hidden', ->
    $(this).prev('span').find('i').toggleClass('icon-chevron-right').toggleClass('icon-chevron-down')


  #$('.collapsible').click ->
  #  visible = !$(this).next().hasClass('collapsed')
  #  if visible
  #    $(this).next().addClass('collapsed')
  #  else
  #    $(this).next().removeClass('collapsed')
  

# vim: ft=coffee foldmethod=manual
