#include <wayfire/plugin.hpp>
#include <wayfire/output.hpp>
#include <wayfire/core.hpp>
#include <wayfire/view.hpp>
#include <wayfire/util/duration.hpp>
#include <wayfire/workspace-manager.hpp>
#include <wayfire/render-manager.hpp>
#include <wayfire/compositor-view.hpp>
#include <wayfire/output-layout.hpp>
#include <wayfire/touch/touch.hpp>
#include <wayfire/plugins/vswitch.hpp>

#include <cmath>
#include <linux/input.h>
#include <wayfire/signal-definitions.hpp>
#include <wayfire/plugins/common/preview-indication.hpp>

#include "snap_signal.hpp"
#include <wayfire/plugins/common/move-snap-helper.hpp>
#include <wayfire/plugins/common/view-change-viewport-signal.hpp>


class wayfire_move3d : public wf::plugin_interface_t
{
    wf::signal_callback_t move_request, view_destroyed;
    wf::button_callback activate_binding;
    wayfire_view view;

    wf::option_wrapper_t<wf::buttonbinding_t> activate_button{"move3d/activate"};

    bool is_using_touch;
    bool was_client_request;

    struct
    {
        nonstd::observer_ptr<wf::preview_indication_view_t> preview;
        int slot_id = 0;
    } slot;

    wf::wl_timer workspace_switch_timer;

#define MOVE_HELPER view->get_data<wf::move_snap_helper_t>()

  public:
    void init() override
    {
        grab_interface->name = "move3d";
        grab_interface->capabilities =
            wf::CAPABILITY_GRAB_INPUT /*| wf::CAPABILITY_MANAGE_DESKTOP*/;

        activate_binding = [=] (uint32_t, int, int)
        {
            is_using_touch     = false;
            was_client_request = false;
            auto view = wf::get_core().get_cursor_focus_view();

            if (view && (view->role != wf::VIEW_ROLE_DESKTOP_ENVIRONMENT))
            {
                return initiate(view);
            }

            return false;
        };

        output->add_button(activate_button, &activate_binding);

        using namespace std::placeholders;
        grab_interface->callbacks.pointer.button =
            [=] (uint32_t b, uint32_t state)
        {
            /* the request usually comes with the left button ... */
            if ((state == WLR_BUTTON_RELEASED) && was_client_request &&
                (b == BTN_LEFT))
            {
                return input_pressed(state, false);
            }

            if (b != wf::buttonbinding_t(activate_button).get_button() /* &&
                b != wf::buttonbinding_t(activate_in3d_button).get_button()*/)
            {
                return;
            }

            is_using_touch = false;
            input_pressed(state, false);
        };

        grab_interface->callbacks.pointer.motion = [=] (int x, int y)
        {
            handle_input_motion();
        };

        grab_interface->callbacks.touch.motion =
            [=] (int32_t id, int32_t sx, int32_t sy)
        {
            handle_input_motion();
        };

        grab_interface->callbacks.touch.up = [=] (int32_t id)
        {
            if (wf::get_core().get_touch_state().fingers.empty())
            {
                input_pressed(WLR_BUTTON_RELEASED, false);
            }
        };

        grab_interface->callbacks.cancel = [=] ()
        {
            input_pressed(WLR_BUTTON_RELEASED, false);
        };

        //move_request =
        //    std::bind(std::mem_fn(&wayfire_move::move_requested), this, _1);
        //output->connect_signal("view-move-request", &move_request);

        view_destroyed = [=] (wf::signal_data_t *data)
        {
            if (get_signaled_view(data) == view)
            {
                input_pressed(WLR_BUTTON_RELEASED, true);
            }
        };
        output->connect_signal("view-disappeared", &view_destroyed);
        output->connect_signal("view-move-check", &on_view_check_move);
    }

    void move_requested(wf::signal_data_t *data)
    {
        auto view = get_signaled_view(data);
        if (!view)
        {
            return;
        }

        auto touch = wf::get_core().get_touch_state();
        is_using_touch = !touch.fingers.empty();

        was_client_request = true;
        initiate(view);
    }

    wf::signal_connection_t on_view_check_move = [=] (wf::signal_data_t *data)
    {
        auto ev = static_cast<wf::view_move_check_signal*>(data);
        if (!ev->can_continue && can_move_view(ev->view))
        {
            ev->can_continue = true;
        }
    };

    /**
     * Calculate plugin activation flags for the view.
     *
     * Activation flags ignore input inhibitors if the view is in the desktop
     * widget layer (i.e OSKs)
     */
    uint32_t get_act_flags(wayfire_view view)
    {
        uint32_t view_layer = output->workspace->get_view_layer(view);
        /* Allow moving an on-screen keyboard while screen is locked */
        bool ignore_inhibit = view_layer == wf::LAYER_DESKTOP_WIDGET;
        uint32_t act_flags  = 0;
        if (ignore_inhibit)
        {
            act_flags |= wf::PLUGIN_ACTIVATION_IGNORE_INHIBIT;
        }

        return act_flags;
    }

    /**
     * Calculate the view which is the actual target of this move operation.
     *
     * Usually, this is the view itself or its topmost parent if the join_views
     * option is set.
     */
    wayfire_view get_target_view(wayfire_view view)
    {
        while (view && view->parent/* && join_views*/)
        {
            view = view->parent;
        }

        return view;
    }

    bool can_move_view(wayfire_view view)
    {
        if (!view || !view->is_mapped())
        {
            return false;
        }

        view = get_target_view(view);

        auto current_ws_impl =
            output->workspace->get_workspace_implementation();
        if (!current_ws_impl->view_movable(view))
        {
            return false;
        }

        return output->can_activate_plugin(grab_interface, get_act_flags(view));
    }

    bool initiate(wayfire_view view)
    {
        view = get_target_view(view);
        if (!can_move_view(view) || (view && (view->get_output() != output)))
        {
            return false;
        }

        if (!output->activate_plugin(grab_interface, get_act_flags(view)))
        {
            return false;
        }

        if (!grab_interface->grab())
        {
            output->deactivate_plugin(grab_interface);

            return false;
        }

        ensure_move_helper_at(view, get_input_coords());

        output->focus_view(view, true);
        /*if (enable_snap)
        {
            slot.slot_id = 0;
        }*/

        this->view = view;
        output->render->set_redraw_always();
        //update_multi_output();

        return true;
    }

    void deactivate()
    {
        grab_interface->ungrab();
        output->deactivate_plugin(grab_interface);
        output->render->set_redraw_always(false);
    }

    void input_pressed(uint32_t state, bool view_destroyed)
    {
        if (state != WLR_BUTTON_RELEASED)
        {
            return;
        }

        deactivate();

        /* The view was moved to another output or was destroyed,
         * we don't have to do anything more */
        if (view_destroyed)
        {
            view->erase_data<wf::move_snap_helper_t>();
            this->view = nullptr;

            return;
        }

        MOVE_HELPER->handle_input_released();

        /* Delete any mirrors we have left, showing an animation */
        //delete_mirror_views(true);

        /* Don't do snapping, etc for shell views */
        if (view->role == wf::VIEW_ROLE_DESKTOP_ENVIRONMENT)
        {
            view->erase_data<wf::move_snap_helper_t>();
            this->view = nullptr;
            return;
        }

#if 0
        if (enable_snap && (slot.slot_id != 0))
        {
            snap_signal data;
            data.view = view;
            data.slot = (slot_type)slot.slot_id;
            output->emit_signal("view-snap", &data);

            /* Update slot, will hide the preview as well */
            update_slot(0);
        }
#endif
        view_change_viewport_signal workspace_may_changed;
        workspace_may_changed.view = this->view;
        workspace_may_changed.to   = output->workspace->get_current_workspace();
        workspace_may_changed.old_viewport_invalid = false;
        output->emit_signal("view-change-viewport", &workspace_may_changed);

        view->erase_data<wf::move_snap_helper_t>();
        this->view = nullptr;
    }

#if 0
    /* Calculate the slot to which the view would be snapped if the input
     * is released at output-local coordinates (x, y) */
    int calc_slot(int x, int y)
    {
        auto g = output->workspace->get_workarea();
        if (!(output->get_relative_geometry() & wf::point_t{x, y}))
        {
            return 0;
        }

        if (view && (output->workspace->get_view_layer(view) != wf::LAYER_WORKSPACE))
        {
            return 0;
        }

        int threshold = snap_threshold;

        bool is_left   = x - g.x <= threshold;
        bool is_right  = g.x + g.width - x <= threshold;
        bool is_top    = y - g.y < threshold;
        bool is_bottom = g.x + g.height - y < threshold;

        bool is_far_left   = x - g.x <= quarter_snap_threshold;
        bool is_far_right  = g.x + g.width - x <= quarter_snap_threshold;
        bool is_far_top    = y - g.y < quarter_snap_threshold;
        bool is_far_bottom = g.x + g.height - y < quarter_snap_threshold;

        int slot = 0;
        if ((is_left && is_far_top) || (is_far_left && is_top))
        {
            slot = SLOT_TL;
        } else if ((is_right && is_far_top) || (is_far_right && is_top))
        {
            slot = SLOT_TR;
        } else if ((is_right && is_far_bottom) || (is_far_right && is_bottom))
        {
            slot = SLOT_BR;
        } else if ((is_left && is_far_bottom) || (is_far_left && is_bottom))
        {
            slot = SLOT_BL;
        } else if (is_right)
        {
            slot = SLOT_RIGHT;
        } else if (is_left)
        {
            slot = SLOT_LEFT;
        } else if (is_top)
        {
            // Maximize when dragging to the top
            slot = SLOT_CENTER;
        } else if (is_bottom)
        {
            slot = SLOT_BOTTOM;
        }

        return slot;
    }

    void update_workspace_switch_timeout(int slot_id)
    {
        if ((workspace_switch_after == -1) || (slot_id == 0))
        {
            workspace_switch_timer.disconnect();

            return;
        }

        int dx = 0, dy = 0;
        if (slot_id >= 7)
        {
            dy = -1;
        }

        if (slot_id <= 3)
        {
            dy = 1;
        }

        if (slot_id % 3 == 1)
        {
            dx = -1;
        }

        if (slot_id % 3 == 0)
        {
            dx = 1;
        }

        if ((dx == 0) && (dy == 0))
        {
            workspace_switch_timer.disconnect();

            return;
        }

        wf::point_t cws = output->workspace->get_current_workspace();
        wf::point_t tws = {cws.x + dx, cws.y + dy};
        wf::dimensions_t ws_dim = output->workspace->get_workspace_grid_size();
        wf::geometry_t possible = {
            0, 0, ws_dim.width, ws_dim.height
        };

        /* Outside of workspace grid */
        if (!(possible & tws))
        {
            workspace_switch_timer.disconnect();

            return;
        }

        workspace_switch_timer.set_timeout(workspace_switch_after, [this, tws] ()
        {
            output->workspace->request_workspace(tws, {this->view});
        });
    }

    void update_slot(int new_slot_id)
    {
        /* No changes in the slot, just return */
        if (slot.slot_id == new_slot_id)
        {
            return;
        }

        /* Destroy previous preview */
        if (slot.preview)
        {
            auto input = get_input_coords();
            slot.preview->set_target_geometry(
                {input.x, input.y, 1, 1}, 0, true);
            slot.preview = nullptr;
        }

        slot.slot_id = new_slot_id;

        /* Show a preview overlay */
        if (new_slot_id)
        {
            snap_query_signal query;
            query.slot = (slot_type)new_slot_id;
            query.out_geometry = {0, 0, -1, -1};
            output->emit_signal("query-snap-geometry", &query);

            /* Unknown slot geometry, can't show a preview */
            if ((query.out_geometry.width <= 0) || (query.out_geometry.height <= 0))
            {
                return;
            }

            auto input   = get_input_coords();
            auto preview = new wf::preview_indication_view_t(output,
                {input.x, input.y, 1, 1});

            wf::get_core().add_view(
                std::unique_ptr<wf::view_interface_t>(preview));

            preview->set_output(output);
            preview->set_target_geometry(query.out_geometry, 1);
            slot.preview = nonstd::make_observer(preview);
        }

        update_workspace_switch_timeout(new_slot_id);
    }
#endif
    /* Returns the currently used input coordinates in global compositor space */
    wf::point_t get_global_input_coords()
    {
        wf::pointf_t input;
        if (is_using_touch)
        {
            auto center = wf::get_core().get_touch_state().get_center().current;
            input = {center.x, center.y};
        } else
        {
            input = wf::get_core().get_cursor_position();
        }

        return {(int)input.x, (int)input.y};
    }

    /* Returns the currently used input coordinates in output-local space */
    wf::point_t get_input_coords()
    {
        auto og     = output->get_layout_geometry();
        auto coords = get_global_input_coords() - wf::point_t{og.x, og.y};

        /* If the currently moved view is not a toplevel view, but a child
         * view, do not move it outside its outpup */
        if (view && view->parent)
        {
            double x   = coords.x;
            double y   = coords.y;
            auto local = output->get_relative_geometry();
            wlr_box_closest_point(&local, x, y, &x, &y);
            coords = {(int)x, (int)y};
        }

        return coords;
    }

#if 0
    struct wf_move_output_state : public wf::custom_data_t
    {
        nonstd::observer_ptr<wf_move_mirror_view> view;
    };

    std::string get_data_name()
    {
        return "wf-move-" + output->to_string();
    }

    /* Delete the mirror view on the given output.
     * If the view hasn't been unmapped yet, then do so. */
    void delete_mirror_view_from_output(wf::output_t *wo,
        bool show_animation, bool already_unmapped)
    {
        if (!wo->has_data(get_data_name()))
        {
            return;
        }

        auto view = wo->get_data<wf_move_output_state>(get_data_name())->view;
        /* We erase so early so that in case of already_unmapped == false,
         * we don't do this again for the unmap signal which will be triggered
         * by our view->unmap() call */
        wo->erase_data(get_data_name());

        view->show_animation = show_animation;
        if (!already_unmapped)
        {
            view->close();
        }

        wo->erase_data(get_data_name());
    }

    /* Destroys all mirror views created by this plugin */
    void delete_mirror_views(bool show_animation)
    {
        for (auto& wo : wf::get_core().output_layout->get_outputs())
        {
            delete_mirror_view_from_output(wo,
                show_animation, false);
        }
    }

    wf::signal_connection_t handle_mirror_view_unmapped =
        [=] (wf::signal_data_t *data)
    {
        auto view = get_signaled_view(data);
        delete_mirror_view_from_output(view->get_output(), true, true);
        view->disconnect_signal(&handle_mirror_view_unmapped);
    };

    /* Creates a new mirror view on output wo if it doesn't exist already */
    void ensure_mirror_view(wf::output_t *wo)
    {
        if (wo->has_data(get_data_name()))
        {
            return;
        }

        auto base_output   = output->get_layout_geometry();
        auto mirror_output = wo->get_layout_geometry();

        auto mirror = new wf_move_mirror_view(view, wo,
            base_output.x - mirror_output.x,
            base_output.y - mirror_output.y);

        wf::get_core().add_view(
            std::unique_ptr<wf::view_interface_t>(mirror));

        auto wo_state = wo->get_data_safe<wf_move_output_state>(get_data_name());
        wo_state->view = nonstd::make_observer(mirror);
        mirror->connect_signal("unmapped", &handle_mirror_view_unmapped);
    }

    /* Update the view position, with respect to the multi-output configuration
     *
     * Views in wayfire are visible on only a single output. However, when the user
     * moves the view between outputs, it is desirable to temporarily show the view
     * on all outputs whose boundaries it crosses. We emulate this behavior by
     * creating
     * mirror views of the view being moved, while fading them in and out when needed
     * */
    void update_multi_output()
    {
        /* We are not in the join_view mode, so we can move dialogues
         * independently of their main view. However, we do not support
         * moving dialogues to a different output than their main view. */
        if (this->view && this->view->parent)
        {
            return;
        }

        /* The mouse isn't on our output anymore -> transfer ownership of
         * the move operation to the other output where the input currently is */
        auto global = get_global_input_coords();
        auto target_output =
            wf::get_core().output_layout->get_output_at(global.x, global.y);

        if (target_output != output)
        {
            /* The move plugin on the next output will create new mirror views */
            delete_mirror_views(false);

            if (wf::can_start_move_on_output(view, target_output))
            {
                /** First, reset moving view so that we don't remove its snap
                 * helper when the output is changed. */
                deactivate();
                auto view_copy = this->view;
                this->view = nullptr;
                wf::start_move_on_output(view_copy, target_output);
            } else
            {
                input_pressed(WLR_BUTTON_RELEASED, false);
            }

            return;
        }

        auto current_og = output->get_layout_geometry();
        auto current_geometry =
            view->get_bounding_box() + wf::point_t{current_og.x, current_og.y};

        for (auto& wo : wf::get_core().output_layout->get_outputs())
        {
            if (wo == output) // skip the same output
            {
                continue;
            }

            auto og = output->get_layout_geometry();
            /* A view is visible on the other output as well */
            if (og & current_geometry)
            {
                ensure_mirror_view(wo);
            }
        }
    }
#endif

    void handle_input_motion()
    {
        auto input = get_input_coords();
        MOVE_HELPER->handle_motion(get_input_coords(), true);

        //update_multi_output();
        /* View might get destroyed when updating multi-output */
#if 0
        if (view)
        {
            // Make sure that fullscreen views are not tiled.
            // We allow movement of fullscreen views but they should always
            // retain their fullscreen state (but they can be moved to other
            // workspaces). Unsetting the fullscreen state can break some
            // Xwayland games.
            if (enable_snap && !MOVE_HELPER->is_view_fixed() &&
                !this->view->fullscreen)
            {
                update_slot(calc_slot(input.x, input.y));
            }
        } else
        {
            /* View was destroyed, hide slot */
            update_slot(0);
        }
#endif
    }

    void fini() override
    {
        if (grab_interface->is_grabbed())
        {
            input_pressed(WLR_BUTTON_RELEASED, false);
            //delete_mirror_views(false);
        }

        output->rem_binding(&activate_binding);
        //output->rem_binding(&activate_in3d_binding);
        //output->disconnect_signal("view-move-request", &move_request);
        output->disconnect_signal("view-disappeared", &view_destroyed);
    }
};

DECLARE_WAYFIRE_PLUGIN(wayfire_move3d);

