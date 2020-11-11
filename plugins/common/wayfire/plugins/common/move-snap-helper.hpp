#pragma once

#include <wayfire/workspace-manager.hpp>
#include <wayfire/view.hpp>
#include <wayfire/core.hpp>
#include <wayfire/option-wrapper.hpp>
#include <wayfire/plugins/wobbly/wobbly-signal.hpp>
#include <wayfire/util/log.hpp>
#include <wayfire/debug.hpp>
#include <wayfire/view-transform.hpp>

namespace wf
{
/**
 * A class which represents the action of moving a view.
 * It provides several conveniences:
 *
 * 1. Interaction with the wobbly plugin
 * 2. Support for locking tiled views in-place until a certain threshold
 * 3. Ensuring view is grabbed at the correct place
 */
constexpr const char *transformer_move3d = "transformer-move3d";
class move_snap_helper_t : public wf::custom_data_t
{
    wayfire_view view;
    wf::point_t grab;

    /* TODO: figure out where to put these options */
    wf::option_wrapper_t<bool> enable_snap_off{"move/enable_snap_off"};
    wf::option_wrapper_t<int> snap_off_threshold{"move/snap_off_threshold"};
    wf::option_wrapper_t<bool> join_views{"move/join_views"};

    bool view_in_slot; /* Whether the view is held at its original position */
    double px, py; /* Percentage of the view width/height from the grab point
                    *  to its upper-left corner */
    bool is_using_3d;

    std::vector<wayfire_view> enum_views(wayfire_view view)
    {
        if (join_views)
        {
            return view->enumerate_views();
        } else
        {
            return {view};
        }
    }

  public:
    move_snap_helper_t(wayfire_view view, wf::point_t grab)
    {
        this->view = view;
        this->grab = grab;
        this->last_grabbing_position = grab;
        this->prev_grabbing_position = grab;

        view_in_slot = should_enable_snap_off();
        for (auto v : enum_views(view))
        {
            start_wobbly(v, grab.x, grab.y);
        }

        auto wmg = view->get_wm_geometry();
        px = 1.0 * (grab.x - wmg.x) / wmg.width;
        py = 1.0 * (grab.y - wmg.y) / wmg.height;
        view->set_moving(1);
        view->connect_signal("geometry-changed", &view_geometry_changed);
    }

    /**
     * Destroy the move snap helper.
     * NB: The destructor will not release the wobbly grab if input has not been
     *   released yet! This is useful if the wobbly grab is to be "transferred"
     *   to another plugin.
     */
    virtual ~move_snap_helper_t()
    {
        view->set_moving(false);
        //view->set_resizing(false);
        view->disconnect_signal("geometry-changed", &view_geometry_changed);
        this->view = nullptr;
    }

    /**
     * Handle a new input event which is a motion of the grabbing point, for ex.
     * cursor or touch point was moved.
     *
     * @param to The new grab point position, in output-local coordinates.
     */
    virtual void handle_motion(wf::point_t to, bool is_using_3d = false)
    {
        this->is_using_3d = is_using_3d;

        if (is_using_3d && !view->get_transformer(transformer_move3d))
        {
            view->add_transformer(std::make_unique<wf::view_3D>(view),
                transformer_move3d);
        }

        for (auto v : enum_views(view))
        {
            move_wobbly(v, to.x, to.y);
        }

        double distance = std::sqrt((to.x - grab.x) * (to.x - grab.x) +
            (to.y - grab.y) * (to.y - grab.y));

        /* Reached threshold */
        if (view_in_slot && (distance >= snap_off_threshold))
        {
            snap_off();
        }

        /* View is stuck, we shouldn't change its geometry */
        if (view_in_slot)
        {
            return;
        }

        this->last_grabbing_position = to;
        adjust_around_grab();
    }

    /**
     * Handle a new input event which is a motion of the grabbing point, for ex.
     * cursor or touch point was moved.
     *
     * In contrast to handle_motion, this event causes wobbly to directly "jump"
     * to the new position.
     *
     * @param to The new grab point position, in output-local coordinates.
     */
    virtual void handle_grab_translated(wf::point_t to)
    {
        wf::point_t delta = to - this->last_grabbing_position;
        for (auto v : enum_views(view))
        {
            translate_wobbly(v, delta);
        }

        this->last_grabbing_position = to;
        adjust_around_grab();
    }

    /**
     * The input point was released (mouse unclicked, finger lifted, etc.)
     * This will also release the wobbly grab.
     */
    virtual void handle_input_released()
    {
        for (auto v : enum_views(view))
        {
            end_wobbly(v);
        }

        //view->pop_transformer(transformer_move3d);
        view->disconnect_signal("geometry-changed", &view_geometry_changed);
        /**
         * Restore the fullscreen state so the window is
         * not in a weird position afterwards, fullscreen
         * are by spec the geometry of the output (or span outputs)
         */
        if (view->fullscreen)
        {
            wf::output_t *output = view->get_output();
            wf::geometry_t output_geometry = output->get_relative_geometry();
            wf::point_t current_ws = output->workspace->get_current_workspace();
            wf::point_t target_ws;
            target_ws.x =
                floor(
                    (last_grabbing_position.x / output_geometry.width) +
                    current_ws.x);
            target_ws.y =
                floor(
                    (last_grabbing_position.y / output_geometry.height) +
                    current_ws.y);
            view->fullscreen_request(output, true, target_ws);
        }
    }

    /** @return Whether the view is freely moving or stays at the same place */
    virtual bool is_view_fixed() const
    {
        return this->view_in_slot;
    }

  protected:
    virtual bool should_enable_snap_off() const
    {
        return enable_snap_off &&
               (view->tiled_edges || view->fullscreen);
    }

    /** Move the view out of its slot */
    virtual void snap_off()
    {
        this->view_in_slot = false;
        if (view->tiled_edges && !view->fullscreen)
        {
            view->tile_request(0);
        }
    }

    wf::point_t prev_grabbing_position;
    wf::point_t last_grabbing_position;
    /** Adjust the view position so that it stays around the grabbing point */
    virtual void adjust_around_grab()
    {
        auto wmg = view->get_wm_geometry();
        wf::point_t target_position = {
            int(last_grabbing_position.x - px * wmg.width),
            int(last_grabbing_position.y - py * wmg.height),
        };

        view->disconnect_signal("geometry-changed", &view_geometry_changed);

        auto transform = dynamic_cast<wf::view_3D*>(
            view->get_transformer(transformer_move3d).get());
        if (is_using_3d && transform)
        {
            glm::vec4 translation_vec = transform->translation_3d * glm::vec4(0.0,0.0,0.0,1.0);
            double translation_z = translation_vec.z;
            double translation_dz = (double)(this->last_grabbing_position.y -
                    this->prev_grabbing_position.y) * -0.5;
            translation_z -= translation_dz;
            if (translation_z > 0.0) translation_z = 0.0;
            if (translation_z < -50.0) translation_z = -50.0;
            view->damage();
            //transform->translation = glm::translate(glm::mat4(1.0), {0.0, 0.0, 0.0});
            transform->translation_3d = glm::translate(glm::mat4(1.0), {0.0, 0.0, translation_z});
            //transform->scaling = glm::scale(glm::mat4(1.0), {1.0f, 1.0f, 1.0});
            //transform->rotation = glm::rotate(glm::mat4(1.0), 0.0f, {0.0, 1.0, 0.0});
            //uint32_t edges = 0;
            //edges |= WLR_EDGE_LEFT;
            //edges |= WLR_EDGE_RIGHT;
            //edges |= WLR_EDGE_TOP;
            //edges |= WLR_EDGE_BOTTOM;
            //view->set_moving(true);
            //view->set_resizing(true, edges);
            view->damage();
        }
        else
        {
            view->move(target_position.x, target_position.y);
        }

        this->prev_grabbing_position = this->last_grabbing_position;
        view->connect_signal("geometry-changed", &view_geometry_changed);
    }

    signal_callback_t view_geometry_changed = [=] (signal_data_t*)
    {
        adjust_around_grab();
    };
};

/**
 * Add a move helper to the view if not already present, otherwise, move the
 * grab to the given point.
 */
void ensure_move_helper_at(wayfire_view view, wf::point_t point)
{
    auto helper = view->get_data<wf::move_snap_helper_t>();
    if (helper != nullptr)
    {
        helper->handle_grab_translated(point);
    } else
    {
        view->store_data(std::make_unique<wf::move_snap_helper_t>(view, point));
    }
}

/**
 * name: view-move-check
 * on: output
 * when: A plugin can emit this signal on an output to check whether there is a
 *   plugin on that output which can continue an interactive move operation.
 */
struct view_move_check_signal : public _view_signal
{
    /** A plugin should set this to true if it can continue a move operation. */
    bool can_continue = false;
};

/**
 * Check whether we can start an interactive move on the other output.
 */
bool can_start_move_on_output(wayfire_view view, wf::output_t *output)
{
    view_move_check_signal check;
    check.view = view;
    output->emit_signal("view-move-check", &check);

    return check.can_continue;
}

/**
 * Start an interactive move on another output.
 * Precondition: the view is being moved with the snap helper.
 */
void start_move_on_output(wayfire_view view, wf::output_t *output)
{
    wf::get_core().move_view_to_output(view, output, false);
    wf::get_core().focus_output(output);
    view->move_request();
}
}
